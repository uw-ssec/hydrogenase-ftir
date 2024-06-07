# I/O Functions for wrapping the ingestion/writing of data.

from brukeropusreader import read_file
import pathlib
from pathlib import Path

def import_run_data(path_to_data:Path):
    """
    Function to import raw spectra data from the file system/a given folder.

    Parameters:
    - path_to_data: Path
        pathlib object which contains the location of the desired data to import.

    Returns:
    - raw_data: dict of OPUS Files
        Dictionary containing the raw data, read in as OPUSFiles

    """
    #Pulling in all sample data from designated file path
    raw_files = list(path_to_data.iterdir())
    raw_files.sort()

    #Initializing dict of raw spectra files from the file system
    raw_data = dict()

    #Populating the raw_test_data dict with all the read in raw opus files
    for i in raw_files:
        if not i.name.startswith('.DS_Store'):
            #FUTURE: pass in results of read_file into a ProSpecPy init function. Then our resulting dict is made up of 
            #ProSpecPy objects, calling the "set_raw"
            #Also passing in output path
            raw_data[f'{i.name}'] = read_file(i)

    print(raw_data.keys())
    return raw_data

