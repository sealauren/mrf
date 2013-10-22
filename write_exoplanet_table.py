#write_exoplanet_table.py
'''Queries exoplanets.org for planets with measured masses and downloads
their parameters.  Adds planets from our research.  Outputs joined table
with all the parameters I want.'''

def query_exoplanets():
    '''Will eventually query exoplanets.org for exoplanets.  Now it just reads in a local
    csv file.'''
    # This is going to require javascript to read the table!  Get help on this.
    import numpy as np
    import math
    import pandas as pd
    import matplotlib.pyplot as plt
    dtype_dict = {'NAME': str, 'MSINI': np.float32, 'MSINIUPPER': np.float32, 'MSINILOWER': np.float32,
                'UMSINI':np.float32, 'R':np.float32, 'RUPPER':np.float32, 'RLOWER':np.float32, 'UR':np.float32,
	      'MSTAR':np.float32, 'MSTARUPPER':np.float32, 'MSTARLOWER':np.float32, 'UMSTAR':np.float32, 
                'RSTAR':np.float32, 'RSTARUPPER':np.float32, 'RSTARLOWER':np.float32, 'URSTAR':np.float32,
              'TEFF':np.float32, 'TEFFUPPER':np.float32, 'TEFFLOWER':np.float32, 'UTEFF':np.float32,
              'A':np.float32, 'AUPPER':np.float32, 'ALOWER':np.float32, 'UA':np.float32,
              'PER':np.float32, 'ECC':np.float32, 'I':np.float32,
              'FIRSTREF':str, 'FIRSTURL':str, 'ORBREF':str, 'ORBURL':str, 'DATE':int}
    df = pd.read_csv("all_exos_errors5.csv",skiprows=[1])
    df.sort(inplace=True)
    
def update_values():
    '''Overwrites exoplanets.org data with more recent, hard-coded data.'''
    pass
    
def add_k22():
    '''Adds planets from the 22 stars observed at Keck/HIRES over the last few
    years.'''
    pass
    
def write_table():
    '''Writes exoplanet properties to a new table.'''
    pass
    
