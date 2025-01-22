import os
import time
import requests
import urllib
import pandas as pd
from tqdm import tqdm
from io import StringIO
from neddy import _basesearch
from ned_extinction_calc import request_extinctions
import json

class PhotometrySearch(_basesearch):
    """
    Perform a NED name-search and return metadata or photometric data for the matched sources.
    """

    BASE_URL = 'https://ned.ipac.caltech.edu/cgi-bin/datasearch?'
    TEMP_FILE_PATH = "/Users/jhzhe/dev/MPhys-Project/output_data/temp/dict.txt"

    def __init__(self, log, names=None, filters=None, ra_list = None, dec_list = None, quiet=False, 
                 verbose=False, output_file_path=None, SED=False):
        """
        Initialize the NameSearch class.

        Args:
            log: Logger instance for logging.
            names: List of names or a single name to search.
            quiet: Suppress stdout output.
            verbose: Provide detailed output.
            output_file_path: Path to save results.
            SED: If True, perform SED-based queries.
        """
        self.log = log
        self.names = [names] if isinstance(names, str) else names
        self.filters = filters
        self.quiet = quiet
        self.ra_list = ra_list
        self.dec_list = dec_list
        self.verbose = verbose
        self.output_file_path = output_file_path
        self.SED = SED

        # Frequency ranges for filtering
        self.FREQ_LOWER = 430 * 10**12  # Hz (red light)
        self.FREQ_UPPER = 800 * 10**12  # Hz (near blue light)
        self.FREQ_MID = 530 * 10**12  # Hz (green light)

        c = 3e8

        # Bandpass dictionary with central wavelengths in micrometers (µm)
        band_dict = {
            "Landolt U": 0.3508,
            "Landolt B": 0.4329,
            "Landolt V": 0.5422,
            "Landolt R": 0.6428
        }

        # Calculate central frequencies (in Hz)
        central_frequencies = {band: c / (wavelength * 1e-6) for band, wavelength in band_dict.items()}

        # Sort bands by central frequency (just in case they're not ordered)
        sorted_bands = sorted(central_frequencies.items(), key=lambda x: x[1], reverse=True)

        # Create continuous frequency ranges
        self.frequency_ranges = {}
        for i in range(len(sorted_bands)):
            band, freq = sorted_bands[i]
            if i == 0:  # First band (highest frequency)
                # Use midpoint between this and the next band as the lower bound
                lower_bound = freq - (freq - sorted_bands[i + 1][1]) / 2
                upper_bound = freq
            elif i == len(sorted_bands) - 1:  # Last band (lowest frequency)
                lower_bound = freq
                upper_bound = freq + (sorted_bands[i - 1][1] - freq) / 2
            else:  # Middle bands
                lower_bound = freq - (freq - sorted_bands[i + 1][1]) / 2
                upper_bound = freq + (sorted_bands[i - 1][1] - freq) / 2
            
            # Save the range for the band
            self.frequency_ranges[band] = (lower_bound, upper_bound)

        # Ensure the temp directory exists
        os.makedirs(os.path.dirname(self.TEMP_FILE_PATH), exist_ok=True)

    def get(self):
        """
        Perform the name search and retrieve results.
        """
        if not self.SED:
            self.log.warning("SED queries are disabled. No data will be fetched.")
            return None

        self.log.info("Starting SED queries for names.")
        return self.perform_sed_queries(self.names)

    def build_sed_query_url(self, name):
        """
        Build the SED query URL for a given name.

        Args:
            name: Name of the object.

        Returns:
            URL for the SED query.
        """
        params = {
            'search_type': 'Photometry',
            'meas_type': 'bot',
            'ebars_spec': 'ebars',
            'label_spec': 'no',
            'x_spec': 'freq',
            'y_spec': 'Fnu_jy',
            'xr': '-1',
            'of': 'ascii_bar',
            'objname': name
        }
        return self.BASE_URL + urllib.parse.urlencode(params)

    def read_processed_names(self):
        """
        Read previously processed names from the temp file.

        Returns:
            List of processed names.
        """
        if not os.path.exists(self.TEMP_FILE_PATH):
            return []

        with open(self.TEMP_FILE_PATH, 'r') as file:
            return file.read().splitlines()

    def save_processed_name(self, name):
        """
        Save a processed name to the temp file.

        Args:
            name: Name to save.
        """
        with open(self.TEMP_FILE_PATH, 'a') as file:
            file.write(name + '\n')
    
    def save_colour_map(self, colour_map):
        """
        Save the colour_map dictionary to the temp file.

        Args:
            colour_map: Dictionary mapping object names to color groups.
        """
        with open(self.TEMP_FILE_PATH, 'w') as file:
            json.dump(colour_map, file)

    def load_colour_map(self):
        """
        Load the colour_map dictionary from the temp file.

        Returns:
            Dictionary containing the colour_map. Returns an empty dictionary if the file does not exist.
        """
        if not os.path.exists(self.TEMP_FILE_PATH):
            return {}

        with open(self.TEMP_FILE_PATH, 'r') as file:
            try:
                return json.load(file)
            except json.JSONDecodeError:
                self.log.warning("Failed to decode JSON. Starting with an empty dictionary.")
                return {}
                
    def download_sed_data(self, name, query_url):
        """
        Download SED data for a single name.

        Args:
            name: Object name.
            query_url: URL for the SED query.

        Returns:
            DataFrame of photometry data or None if no valid data is found.
        """
        try:
            response = requests.get(query_url, timeout=10)
            response.raise_for_status()

            # Parse the response into a DataFrame
            data = pd.read_csv(StringIO(response.content.decode('utf-8')),
                               sep='|',
                               skiprows=16,
                               engine='python')

            return data
        except requests.exceptions.RequestException as e:
            self.log.error(f"Network error for {name}: {e}")
            return None
        except pd.errors.EmptyDataError:
            self.log.warning(f"No valid data found for {name}.")
            return None
        except Exception as e:
            self.log.error(f"Unexpected error for {name}: {e}")
            return None

    def process_photometry_data(self, df, ra, dec):
        """
        Process photometry data to determine the object's "color group."

        Args:
            name: Object name.
            df: DataFrame containing photometry data.
            ra: Right ascension of the object. List
            dec: Declination of the object. List

        Returns:
            Color group for the object ('1', '2', or '3').
        """
        # Filter by frequency range
        df = df[(df['Frequency'] >= self.FREQ_LOWER) & (df['Frequency'] <= self.FREQ_UPPER)]

        # Extinction values
        extinctions = request_extinctions(ra, dec, filters=self.filters, as_dict=True)
        # Example DataFrame filtering based on frequency ranges
        for band, (low_freq, high_freq) in self.frequency_ranges.items():
            df = df.copy()
            df.loc[(df['Frequency'] >= low_freq) & (df['Frequency'] <= high_freq), 'Band'] = band

        if extinctions:
            for key, value in extinctions.items():
                if key in df['Band'].values:
                    df.loc[df['Band'] == key, 'Flux Density Adjusted'] = df['Flux Density'] * 10**(-value / 2.5)

        df['Flux Density Adjusted'] = df['Flux Density Adjusted'].fillna(df['Flux Density'])

        if not df.empty:
            max_flux_freq = df.nlargest(1, 'Flux Density Adjusted')['Frequency'].values[0]
            if self.FREQ_UPPER - 200 * 10**12 <= max_flux_freq <= self.FREQ_UPPER:
                return '1'
            elif self.FREQ_MID <= max_flux_freq < self.FREQ_UPPER - 200 * 10**12:
                return '2'
        return '3'

    def perform_sed_queries(self, names):
        """
        Perform SED queries for a list of names.

        Args:
            names: List of object names.

        Returns:
            Dictionary mapping names to their "color group."
        """
        colour_map = self.load_colour_map()
        processed_names = set(colour_map.keys())

        for name in tqdm(names, desc="Downloading SED data"):
            if name in processed_names:
                self.log.info(f"{name} already processed. Skipping...")
                continue
            
            query_url = self.build_sed_query_url(name)
            ra = self.ra_list[names.index(name)]
            dec = self.dec_list[names.index(name)]
            df = self.download_sed_data(name, query_url)
            time.sleep(1)
            if df is not None:
                colour_map[name] = self.process_photometry_data(df, ra, dec)
            else:
                colour_map[name] = 'No Data'

            self.save_colour_map(colour_map)
            time.sleep(1)  # Avoid overwhelming the server

        return colour_map
