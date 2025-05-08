from __future__ import print_function

import numpy as np
from astroquery.linelists.cdms import CDMS
import astropy.units as u
import pandas as pd
from scipy import constants

pd.options.mode.chained_assignment = None

def generate_einstein_As(elo,lgint,freq,gup,Q):
    """
    Calculate Einstein A for a given transition in a CDMS database
    Follows the calculation detailed in 
    https://cdms.astro.uni-koeln.de/classic/predictions/description.html
    """
    import astropy.units as u
    import astropy.constants as const
    freq = freq*u.MHz

    #Convert elo from cm-1 to K and calculate eup
    E_l = (elo / u.cm) * (const.h * (const.c.to(u.cm/u.s))) / const.k_B
    E_u = E_l + (f.to(1/u.s) * const.h / const.k_B)

    intensity = 10**lgint * u.nm**2 * u.MHz

    gl_S_mu2 = 2.50251e4 * intensity.value * Q['300.0']  / (freq * (np.exp(-E_l / (300*u.K)) - np.exp(-E_u / (300*u.K))))

    EA = 1.16395e-20 * (freq**3) * gl_S_mu2 / gup

    return EA, E_u.value

def getCDMS(mol, emax, fmin=0*u.GHz, fmax=5000*u.GHz, saveFile=True):
    """
    Import CDMS (or JPL) molecular data. 
    Convert it to a LAMDA-style linelist for a given Emax.

    Inputs:
    mol - CDMS (or JPL) molecule index number
    emax - maximum energy level (cm-1) to be included in the linelist
    fmin - lowest transition frequency desired
    fmax - highest transition frequency desired
    saveFile - whether or not to save the linelist to a text file

    Outputs:
    df - dataframe containing molecular parameters, optionally written to file
    """

    #Retrieve the linelist
    ctbl = CDMS.query_lines(min_frequency=fmin, max_frequency=fmax, molecule=mol).to_pandas()

    #Drop NA columns to get rid of empty quantum numbers
    ctbl = ctbl.dropna(axis=1)

    #Drop transitions with Elo > Emax
    ctbl = ctbl[ctbl['ELO'] <= emax]

    #Drop 'name' and 'Lab' columns
    ctbl = ctbl.drop(['name','Lab'], axis=1)

    #Get the molecule data
    molData = CDMS.get_species_table()
    mtbl = molData[molData['tag'] == int(mol.split()[0])]
    Q300 = 10**mtbl['lg(Q(300))'].item()

    ctbl['EA'] = ctbl.apply(lambda x: generate_Einstein_As(x['ELO'],x['LGINT'],x['FREQ'],x['GUP'],Q300), axis=1)

    #Work out how many quantum numbers we have. They come after QNFMT, which is column 8. 
    #Since we have upper and lower, divide by two to figure this out
    cols = ctbl.columns
    nQuanta = int(len(cols[9:]))/2.

    #Generate a list of levels and statistical weights
    enlevels = {}
    gs = {}
    i=0
    for elo, gup, qn2 in zip(ctbl['ELO'],ctbl['GUP']):
        enlevels[]

    

