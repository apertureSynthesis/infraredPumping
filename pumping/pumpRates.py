# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
from astroquery.lamda import Lamda
from astroquery.hitran import Hitran
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
        self.tbl = Hitran.query_lines(molecule_number=id, isotopologue_number=iso, min_frequency=0./u.cm, max_frequency=numax).to_pandas()

        #Translate the bytes to strings
        self.tbl['global_upper_quanta']=self.tbl['global_upper_quanta'].apply(lambda x: x.decode("utf-8"))
        self.tbl['global_lower_quanta']=self.tbl['global_lower_quanta'].apply(lambda x: x.decode("utf-8"))
        self.tbl['local_lower_quanta']=self.tbl['local_lower_quanta'].apply(lambda x: x.decode("utf-8"))
        self.tbl['local_upper_quanta']=self.tbl['local_upper_quanta'].apply(lambda x: x.decode("utf-8"))