# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
from astroquery.lamda import Lamda
from astroquery.hitran import Hitran
from radis.io import fetch_geisa
import astropy.units as u
import pandas as pd
from scipy import constants
import itertools
import os
from infraredPumping.utils.molParams import *
from infraredPumping.utils.lamdaRoutines import *
from infraredPumping.utils.quantumDicts import *

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
        lamdamol, id, iso, includeLevels, groundState, species = getMolParams(self.mol,self.levels)
    
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

        #Import the LAMDA data
        self.collrates, self.radtransitions, self.enlevels = Lamda.query(mol=lamdamol)
        self.levels = self.enlevels.to_pandas()