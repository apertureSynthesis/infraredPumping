# -*- coding: utf-8 -*-

from __future__ import print_function

import os,sys,collections
import numpy as np
from astroquery.hitran import Hitran
from astroquery.linelists.cdms import CDMS
from radis.io import fetch_geisa
import astropy.units as u
import pandas as pd
from scipy import constants
import itertools
from infraredPumping.utils import molParams
from infraredPumping.utils import lamdaRoutines
from infraredPumping.utils import quantumNumbers
from infraredPumping.utils import cdmsRoutines

pd.options.mode.chained_assignment = None
MINTRANS = 10

def nest_defaultdict():
    """ I couldn't get defaultdict to cooperate without making this nested structure.
        Now it works.
    """
    return collections.defaultdict(nest_defaultdict)

class pumpRates(object):
    """
    Calculate effective infrared pumping rates excited by solar radiation as a
    black body radiation field
    """

    def __init__(self, mol, emax, numax = 30000./u.cm, nlev=-1, levels='all'):
        """
        pumping rates excited by solar radiation ignoring hot band cascades
        """

        self.mol = mol
        self.emax = emax
        self.numax = numax
        self.levels = levels

        #Get pertinent information about the molecule
        cdmsMol, id, iso, includeLevels, groundState, species = molParams.getMolParams(self.mol,self.levels)
        print(includeLevels)
    
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
        self.tbl[['local_lower_quanta_clean','local_upper_quanta_clean']] = self.tbl.apply(lambda x: pd.Series(quantumNumbers.generate_linear_quanta(x['local_lower_quanta'],x['local_upper_quanta'],self.mol)), axis=1)

        #Map the HITRAN local quanta to CDMS format for mapping between the two
        self.tbl[['local_lower_quanta_cdms','local_upper_quanta_cdms']] = self.tbl.apply(lambda x: pd.Series(quantumNumbers.hitran_to_cdms(x['local_lower_quanta_clean'],x['local_upper_quanta_clean'],self.mol)), axis=1)

        #Calculate statistical weights for HNC, as they aren't included in GEISA
        if self.mol == 'HNC':
            self.tbl.loc[:,'gp'] = self.tbl['local_uppr_quanta_lamda'].apply(quantumNumbers.generate_statistical_weights)
            self.tbl.loc[:,'gpp'] = self.tbl['local_lower_quanta_lamda'].apply(quantumNumbers.generate_statistical_weights)

        # #Import the LAMDA data
        # collrates, radtransitions, enlevels = lamdaRoutines.query_lamda(lamdamol=lamdamol)
        # self.levels = enlevels.to_pandas()

        #Import the CDMS data
        #First the line list
        self.ctbl, self.enlevels = cdmsRoutines.getCDMS(mol = cdmsMol, emax = self.emax, min_frequency=0*u.GHz, max_frequency=5000*u.GHz)

        #Make a dictionary to hold all of the level rate summations
        self.gratesum={} 

        #Select the upper quantum levels that we want to work with
        vups = self.tbl['global_upper_quanta'].unique().tolist()
        vups.remove(groundState)

        try:
            vups.remove('               ')
        except:
            pass

        if includeLevels != None:
            vups = includeLevels

        pumprates = []
        vuplist = []

        for vup in vups:
            sys.stderr.write( "\nReading transitions for v_up = "+vup.strip()+":\n")
            self.tblf=self.tbl[self.tbl['global_upper_quanta'].isin([vup]) & self.tbl['global_lower_quanta'].isin([groundState])]

            #Arrays for holding the quantum numbers, Einstein A's and g rates
            self.Jinits  = []
            self.Jups    = []
            self.Jfinals = []
            self.As      = []
            self.grates  = []

            #Arrays for wavenumber, statistical weights, and Einstein A
            nu = self.tblf['nu']
            glo= self.tblf['gpp']
            gup= self.tblf['gp']
            A  = self.tblf['a']

            #Get the g factors
            grate = quantumNumbers.glu(nu,gup,glo)

            #Pull out the lower and upper quanta for each transition
            QNL = np.asarray([s for s in self.tblf['local_lower_quanta_cdms']]).transpose() #lower quantum numbers
            QNU = np.asarray([s for s in self.tblf['local_upper_quanta_cdms']]).transpose() #upper quantum numbers

            print('Browsing through the transitions...')
            for Qlo,Qup,A1,G in zip(QNL,QNU,A,grate):
                #For a given Qlo,Qup find other transitions for which Qup is the lower quantum number
                #These represent possible vibronic transitions
                for i in self.tblf[self.tblf['local_upper_quanta_cdms'].str.match(Qup)]['local_lower_quanta_cdms'].values:
                    self.Jfinals.append(i)
                #Record Qlo, Qup, A, and g for every possible vibronic transition
                for j in range(len(self.tblf[self.tblf['local_upper_quanta_cdms'].str.match(Qup)]['local_lower_quanta_cdms'].values)):
                    self.Jinits.append(Qlo)
                    self.Jups.append(Qup)
                    self.As.append(A1)
                    self.grates.append(G)

            if len(self.tblf)>MINTRANS:
                #Remove transitions not in Jup,Jlo
                todel1=np.where(np.asarray([trans in list(zip(QNL,QNU)) for trans in zip(self.Jfinals,self.Jups)])==False)
                todel2=np.where(np.asarray([trans in list(zip(QNL,QNU)) for trans in zip(self.Jinits,self.Jups)])==False)
                self.todel=np.concatenate([todel1[0],todel2[0]])

                self.Jinits  = np.delete(self.Jinits,self.todel)
                self.Jfinals = np.delete(self.Jfinals,self.todel)
                self.Jups    = np.delete(self.Jups,self.todel)
                self.As      = np.delete(self.As,self.todel)
                self.grates  = np.delete(self.grates,self.todel)
                print('Creating the Einstein-A table')
                #Einstein A lookup table
                self.Atable = collections.defaultdict(dict)
                self.Asum = {}   

                for Q1 in self.Jups:
                    self.Asum[Q1] = 0.0
                for Q1,Q2,A12 in zip(QNU,QNL,A):
                    #self.Atable[Q1] = {}
                    self.Atable[Q1][Q2] = A12
                    self.Asum[Q1] += A12


                self.prob = []
                #Fraction of each Jupper going into each Jfinal
                for Q1,Q2 in zip(self.Jups,self.Jfinals):
                    self.prob.append(self.Atable[Q1][Q2] / self.Asum[Q1])

                print(self.enlevels)
                #Multiply by the probabilities and Einstein A's
                self.gratesA = np.array(self.grates) * np.array(self.As) * np.array(self.prob)
                #Sum rates with the same Jinit, Jfinal
                for Jinit,Jfinal,grateJ in zip(self.Jinits,self.Jfinals,self.gratesA):
                    if (Jinit != Jfinal):
                        if (any(self.enlevels.J == Jinit) and any(self.enlevels.J == Jfinal)):
                            try:
                                self.gratesum[str(Jinit+" "+Jfinal)] += grateJ
                            except:
                                self.gratesum[str(Jinit+" "+Jfinal)] = grateJ
            
                vuplist.append(vup)

            else:
                sys.stderr.write("Insufficient transitions to calculate pumping for this level\n")

        # print('grate = ....')
        #print(self.gratesum)

        #Get the unique set Qinit, Qfinal and print the results. Format them so that Qinit, Qfinal correspond to their indices in the LAMDA file
        self.gratesumfinal={}
        QinitSum,QfinalSum = np.asarray([s.split(' ') for s in list(self.gratesum.keys())],dtype=str).transpose()
        sys.stderr.write("\nTotal effective %s pumping rates (Qi, Qf, Rate (s^-1) at 1 AU:\n" %(mol))
        levels_lo = []
        levels_up = []

        for Qi,Qf in sorted(zip(QinitSum,QfinalSum)):
            level_lo = self.enlevels[self.enlevels['J'].str.match(Qi)]['Level'].values[0]
            level_up = self.enlevels[self.enlevels['J'].str.match(Qf)]['Level'].values[0]
            levels_lo.append(level_lo)
            levels_up.append(level_up)
            self.gratesumfinal[str(level_lo)+" "+str(level_up)]=self.gratesum[str(Qi)+" "+str(Qf)]
        file_out = "g_"+str(self.mol)+"_1au.dat"
        f = open(file_out,"w")
        for Li,Lf in sorted(zip(levels_lo,levels_up)):
            f.write('%s %s %e \n' %(Li,Lf,self.gratesumfinal[str(Li)+" "+str(Lf)]))
        f.close()




