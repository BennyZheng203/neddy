#!/usr/local/bin/python
# encoding: utf-8
"""
*Perform a NED name-search and return the metadata for the matched sources*

:Author:
    David Young
"""
from __future__ import print_function
from neddy import _basesearch
import urllib
import os
import sys
from astropy.io.votable import parse_single_table
from io import BytesIO
import pandas as pd
from future import standard_library
standard_library.install_aliases()
import requests
from io import StringIO
import time
import tqdm

class namesearch(_basesearch):
    """
    *Perform a NED name-search and return the metadata for the matched sources*

    **Key Arguments**

    - ``log`` -- logger
    - ``name`` -- name
    - ``quiet`` -- don't print to stdout
    - ``verbose`` -- return more metadata for matches
    - ``searchParams`` -- list of dictionaries to prepend to results
    - ``outputFilePath`` -- path to file to output results to

    **Usage**

    ```python
    from neddy import namesearch
    search = namesearch(
        log=log,
        names=objectName,
        verbose=True,
        outputFilePath="/path/to/output.csv"
    )
    results = search.get()
    ```

    """

    def __init__(
            self,
            log,
            names=False,
            quiet=False,
            verbose=False,
            searchParams=False,
            outputFilePath=False,
            SED=False
    ):
        self.log = log
        log.debug("instantiating a new 'namesearch' object")
        self.names = names
        self.quiet = quiet
        self.verbose = verbose
        self.searchParams = searchParams
        self.outputFilePath = outputFilePath
        self.SED = SED

        self.greenband = 530*10**12

        # CREATE A LIST IF SINGLE NAME GIVEN
        os.environ['TERM'] = 'vt100'
        if not isinstance(self.names, list):
            self.uplist = [self.names]
        else:
            self.uplist = self.names

        return None

    def get(self):
        """
        *perform NED name searches and return the results*

        **Return**

        - ``results`` -- the search results (list of dictionaries)
        """
        self.log.debug('starting the ``get`` method')

        # SPLIT THE LIST OF NAMES INTO BATCHES
        self.theseBatches, self.theseBatchParams = self._split_incoming_queries_into_batches(
            sources=self.uplist,
            searchParams=self.searchParams
        )

        # PERFORM NAME QUERIES AGAINST NED
        if not self.SED:
            self._build_api_url_and_download_results()
            self.results, self.headers = self._parse_the_ned_list_results()
            self._output_results()
        else:
            '''SED query which returns a colour map of the objects wrt to the photometry max flux in visual spectral region'''
            colour_map = self.download_sed_query(self.names)
            self.log.debug('complted the ``get`` method')
            return colour_map

        self.log.debug('completed the ``get`` method')
        return self.results

    def _build_api_url_and_download_results(
            self):
        """
        *build an API URL call for NED to perform batch name queries and download results*
        """
        self.log.debug(
            'completed the ````_build_api_url_and_download_results`` method')

        import urllib.parse
        from fundamentals.download import multiobject_download

        baseUrl = "https://ned.ipac.caltech.edu/cgi-bin/"
        command = "gmd"
        urlParameters = {
            "delimiter": "bar",
            "NO_LINKS": "1",
            "nondb": ["row_count", "user_name_msg", "user_objname"],
            "crosid": "objname",
            "enotes": "objnote",
            "position": ["ra,dec", "bhextin", "pretype", "z", "zunc", "zflag"],
            "gadata": ["magnit", "sizemaj", "sizemin", "morphol"],
            "attdat_CON": ["M", "S", "H", "R", "z"],
            "distance_CON": ["mm", "dmpc"],
            "attdat": "attned"
        }

        queryBase = "%(baseUrl)s%(command)s?uplist=" % locals()
        queryList = []

        # BUILD THE LIST OF QUERIES
        for batch in self.theseBatches:
            thisLength = len(batch)
            queryUrl = queryBase
            # ADD NAMES
            for thisIndex, thisName in enumerate(batch):
                queryUrl = queryUrl + urllib.parse.quote(thisName)
                if thisIndex < thisLength - 1:
                    queryUrl = queryUrl + "%0D"
            # ADD PARAMETERS
            for k, v in list(urlParameters.items()):
                if isinstance(v, list):
                    for item in v:
                        queryUrl = queryUrl + "&" + \
                            k + "=" + urllib.parse.quote(item)
                else:
                    queryUrl = queryUrl + "&" + k + "=" + urllib.parse.quote(v)
            queryList.append(queryUrl)

        # PULL THE RESULT PAGES FROM NED
        self.nedResults = multiobject_download(
            urlList=queryList,
            downloadDirectory="/tmp",
            log=self.log,
            timeStamp=1,
            timeout=3600,
            concurrentDownloads=10,
            resetFilename=False,
            credentials=False,  # { 'username' : "...", "password", "..." }
            longTime=True,
            indexFilenames=True
        )

        for thisIndex, r in enumerate(self.nedResults):
            if r == None:
                thisUrl = queryList[thisIndex]
                self.log.error(
                    'cound not download NED results for URL %(thisUrl)s' % locals())
                sys.exit(0)

        self._convert_html_to_csv()

        self.log.debug(
            'completed the ``_build_api_url_and_download_results`` method')
        return None

    def _output_results(
            self):
        """
        *print the NED search results to STDOUT and/or an output file*
        """
        self.log.debug('starting the ``_output_results`` method')

        from fundamentals.renderer import list_of_dictionaries

        if len(self.results) == 0:
            content = "No results found"
            csvContent = "No results found"
        else:
            dataSet = list_of_dictionaries(
                log=self.log,
                listOfDictionaries=self.results,
            )
            content = dataSet.table(filepath=None)
            if self.outputFilePath:
                csvContent = dataSet.csv(filepath=self.outputFilePath)

        if self.quiet == False:
            print(content)

        self.log.debug('completed the ``_output_results`` method')
        return None
    
    def build_sed_query_url(self, names):
        base_url = 'https://ned.ipac.caltech.edu/cgi-bin/datasearch?'
        params = {
            'search_type': 'Photometry',
            'meas_type': 'bot',
            'ebars_spec': 'ebars',
            'label_spec': 'no',
            'x_spec': 'freq',
            'y_spec': 'Fnu_jy',
            'xr': '-1',
            'of': 'ascii_bar'
        }

        return [base_url + urllib.parse.urlencode({**params, 'objname': name}) for name in names]
    
    def download_sed_query_1(self, names):
        queryList = self.build_sed_query_url(names)

        colour_map = {}
        # Download photometry data for each object, determine colour map and save to file. 
        if not os.path.exists("/Users/jhzhe/dev/MPhys-Project/output_data/temp/dict.txt"):
            f = open("/Users/jhzhe/dev/MPhys-Project/output_data/temp/dict.txt","w+")
            f.close()

        f = open("/Users/jhzhe/dev/MPhys-Project/output_data/temp/dict.txt","r+")
        for name, query in tqdm(zip(names, queryList), total=len(names), desc="Downloading SED data"):
                if name in f.read():
                    f.close()
                    continue
                
                response = requests.get(query)
                data = response.content.decode('utf-8')

                df = pd.read_csv(
                    StringIO(data),
                    sep='|',
                    skiprows=16,
                    engine='python'
                )


                df = df[(df['Frequency'] >= 430 * 10**12) & (df['Frequency'] <= 800 * 10**12)]
                if df.nlargest(1, 'Flux Density')['Flux Density'].values[0] <= 800*10**12 and df.nlargest(1, 'Flux Density')['Flux Density'].values[0] >= 600*10**12:
                    colour_map[name] = '1'
                elif df.nlargest(1, 'Flux Density')['Flux Density'].values[0] >= 530*10**12 and df.nlargest(1, 'Flux Density')['Flux Density'].values[0] <= 600*10**12:
                    colour_map[name] = '2'
                else:
                    colour_map[name]= '3'
                
                f.seek(0)
                f.truncate()
                f.write(str(dict))
                time.sleep(1) # sleep for 1 second to avoid ban
                    
        return colour_map
    
    def download_sed_query(self, names):
        query_list = self.build_sed_query_url(names)

        colour_map = {}
        temp_file_path = "/Users/jhzhe/dev/MPhys-Project/output_data/temp/dict.txt"

        # Ensure the directory exists
        os.makedirs(os.path.dirname(temp_file_path), exist_ok=True)

        try:
            # Ensure the temp file exists
            if not os.path.exists(temp_file_path):
                with open(temp_file_path, "w+") as f:
                    pass

            with open(temp_file_path, "r+") as f:
                processed_names = f.read().splitlines()

                for name, query in tqdm(zip(names, query_list), total=len(names), desc="Downloading SED data"):
                    try:
                        if name in processed_names:
                            continue

                        # Fetch the data
                        response = requests.get(query, timeout=10)
                        response.raise_for_status()  # Raise an HTTPError for bad responses
                        data = response.content.decode('utf-8')

                        # Parse the data
                        df = pd.read_csv(
                            StringIO(data),
                            sep='|',
                            skiprows=16,
                            engine='python'
                        )

                        # Filter frequencies and determine color mapping
                        df = df[(df['Frequency'] >= 430 * 10**12) & (df['Frequency'] <= 800 * 10**12)]
                        
                        if not df.empty:
                            max_flux = df.nlargest(1, 'Flux Density')['Flux Density'].values[0]

                            if 600 * 10**12 <= max_flux <= 800 * 10**12:
                                colour_map[name] = '1'
                            elif 530 * 10**12 <= max_flux < 600 * 10**12:
                                colour_map[name] = '2'
                            else:
                                colour_map[name] = '3'
                        else:
                            colour_map[name] = 'No Data'

                        # Log the processed name
                        f.write(name + '\n')
                        f.flush()

                    except requests.exceptions.RequestException as e:
                        print(f"Network error for {name}: {e}")
                        continue
                    except pd.errors.EmptyDataError:
                        print(f"Data parsing error for {name}: No valid data found.")
                        colour_map[name] = 'No Data'
                        continue
                    except Exception as e:
                        print(f"Unexpected error for {name}: {e}")
                        continue

                    # Avoid overwhelming the server
                    time.sleep(1)

        except IOError as e:
            print(f"File I/O error: {e}")
        except Exception as e:
            print(f"Unexpected error: {e}")

        return colour_map
