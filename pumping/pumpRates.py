# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
from astroquery.hitran import Hitran
from astroquery.linelists.cdms import CDMS
from radis.io import fetch_geisa
import astropy.units as u
import pandas as pd
from scipy import constants
import itertools
import os
from infraredPumping.utils import molParams
from infraredPumping.utils import lamdaRoutines
from infraredPumping.utils import quantumNumbers

pd.options.mode.chained_assignment = None

class pumpRates(object):
    """
    Calculate effective infrared pumping rates excited by solar radiation as a
    black body radiation field
    """

    def __init__(self, mol, numax = 30000./u.cm, nlev=-1, levels='all'):
        """
        pumping rates excited by solar radiation ignoring hot band cascades
        """

        self.mol = mol
        self.numax = numax
        self.levels = levels

        #Get pertinent information about the molecule
        lamdamol, id, iso, includeLevels, groundState, species = molParams.getMolParams(self.mol,self.levels)
    
        #Table of all transitions from HITRAN between the specified min and max frequencies
        #Query GEISA for HNC, otherwise query HITRAN
        if self.mol == 'HNC':
            self.tbl = fetch_geisa(molecule=self.mol,load_wavenum_min=0, load_wavenum_max=self.numax.value)
            #Rename columns to match with HITRAN format
            cnames = ['nu','b','c','elower','global_upper_quanta','global_lower_quanta','local_upper_quanta','local_lower_quanta','f','geisa_iso_id','geisa_mol_id','j','hitran_mol_id','hitran_iso_id','a','n','o','r','ap','bp','cp','fp','op','rp','np','s','sp','t','tp','u','up']
            key_dict = dict([(i,j) for i,j in zip(self.tbl.keys(),cnames)])
            self.tbl.rename(columns=key_dict,inplace=True)
        else:
            self.tbl = Hitran.query_lines(molecule_number=id, isotopologue_number=iso, min_frequency=0./u.cm, max_frequency=numax).to_pandas()

        #Translate the bytes to strings
        self.tbl['global_upper_quanta']=self.tbl['global_upper_quanta'].apply(lambda x: x.decode("utf-8"))
        self.tbl['global_lower_quanta']=self.tbl['global_lower_quanta'].apply(lambda x: x.decode("utf-8"))
        self.tbl['local_lower_quanta']=self.tbl['local_lower_quanta'].apply(lambda x: x.decode("utf-8"))
        self.tbl['local_upper_quanta']=self.tbl['local_upper_quanta'].apply(lambda x: x.decode("utf-8"))

        #Split off the spin species (if applicable) to save time
        self.tbl = quantumNumbers.select_species(self.tbl, self.mol)

        #Generate upper quanta for linear molecules if necessary
        self.tbl['local_lower_quanta_clean','local_upper_quanta_clean'] = self.tbl.apply(lambda x: pd.Series(quantumNumbers.generate_linear_quanta(x['local_lower_quanta'],x['local_upper_quanta'],self.mol), axis=1))

        #Map the HITRAN local quanta to CDMS format for mapping between the two
        self.tbl[['local_lower_quanta_cdms','local_upper_quanta_cdms']] = self.tbl.apply(lambda x: pd.Series(quantumNumbers.hitran_to_cdms(x['local_lower_quanta'],x['local_upper_quanta'],self.mol)), axis=1)

        #Calculate statistical weights for HNC, as they aren't included in GEISA
        if self.mol == 'HNC':
            self.tbl.loc[:,'gp'] = self.tbl['local_uppr_quanta_lamda'].apply(quantumNumbers.generate_statistical_weights)
            self.tbl.loc[:,'gpp'] = self.tbl['local_lower_quanta_lamda'].apply(quantumNumbers.generate_statistical_weights)

        #Import the LAMDA data
        collrates, radtransitions, enlevels = lamdaRoutines.query_lamda(lamdamol=lamdamol)
        self.levels = enlevels.to_pandas()

        #Make a dictionary to hold all of the level rate summations
        gratesum={} 

        #Select the upper quantum levels that we want to work with
        vups = self.tbl['global_upper_quanta'].unique().tolist()

        vups.remove(groundState)

        try:
            vups.remove('               ')
        except:
            pass

        try:
            vups = includeLevels
        except:
            pass

