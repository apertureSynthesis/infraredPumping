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

    def __init__(self, mol, numax = 30000./u.cm, nlev=-1, levels='all'):
        """
        pumping rates excited by solar radiation ignoring hot band cascades
        """

        self.mol = mol
        self.numax = numax
        self.levels = levels

        #Get pertinent information about the molecule
        cdmsMol, id, iso, includeLevels, groundState, species = molParams.getMolParams(self.mol,self.levels)
    
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
        self.ctbl, self.enlevels = cdmsRoutines.getCDMS(mol = cdmsMol, emax = 1000./u.cm, min_frequency=0*u.GHz, max_frequency=5000*u.GHz)

        #Make a dictionary to hold all of the level rate summations
        gratesum={} 

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
            Jinits  = []
            Jups    = []
            Jfinals = []
            As      = []
            grates  = []

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

            for Qlo,Qup,A1,G in zip(QNL,QNU,A,grate):
                #For a given Qlo,Qup find other transitions for which Qup is the lower quantum number
                #These represent possible vibronic transitions
                for i in self.tblf[self.tblf['local_upper_quanta_cdms'].str.match(Qup)]['local_lower_quanta_cdms'].values:
                    Jfinals.append(i)
                #Record Qlo, Qup, A, and g for every possible vibronic transition
                for j in range(len(self.tblf[self.tblf['local_upper_quanta_cdms'].str.match(Qup)]['local_lower_quanta_cdms'].values)):
                    Jinits.append(Qlo)
                    Jups.append(Qup)
                    As.append(A1)
                    grates.append(G)

            if len(self.tblf)>MINTRANS:
                #Remove transitions not in Jup,Jlo
                todel1=np.where(np.asarray([trans in list(zip(QNL,QNU)) for trans in zip(Jfinals,Jups)])==False)
                todel2=np.where(np.asarray([trans in list(zip(QNL,QNU)) for trans in zip(Jinits,Jups)])==False)
                todel=np.concatenate([todel1[0],todel2[0]])

                Jinits  = np.delete(Jinits,todel)
                Jfinals = np.delete(Jfinals,todel)
                Jups    = np.delete(Jups,todel)
                As      = np.delete(As,todel)
                grates  = np.delete(grates,todel)

                #Einstein A lookup table
                Atable = collections.defaultdict(nest_defaultdict)
                Asum = {}         

                for Q1 in Jups:
                    Asum[Q1] = 0.0
                for Q1,Q2,A12 in zip(QNU,QNL,A):
                    try:
                        Atable[Q1][Q2] = A12
                        Asum[Q1] += A12
                    except:
                        pass

                prob = []
                #Fraction of each Jupper going into each Jfinal
                for Q1,Q2 in zip(Jups,Jfinals):
                    prob.append(Atable[Q1][Q2] / Asum[Q1])

                #Multiply by the probabilities and Einstein A's
                gratesA = np.array(grates) * np.array(As) * np.array(prob)
                print(self.enlevels)
                print(Jinits)
                #Sum rates with the same Jinit, Jfinal
                for Jinit,Jfinal,grateJ in zip(Jinits,Jfinals,gratesA):
                    if (Jinit != Jfinal):
                        if (any(self.enlevels.J) == Jinit) and (any(self.enlevels.J == Jfinal)):
                            try:
                                gratesum[str(Jinit+" "+Jfinal)] += grateJ
                            except:
                                gratesum[str(Jinit+" "+Jfinal)] = grateJ
            
                vuplist.append(vup)

            else:
                sys.stderr.write("Insufficient transitions to calculate pumping for this level\n")

        print('grate = ....')
        print(gratesum)

        #Get the unique set Qinit, Qfinal and print the results. Format them so that Qinit, Qfinal correspond to their indices in the LAMDA file
        # gratesumfinal={}
        # QinitSum,QfinalSum = np.asarray([s.split(' ') for s in list(gratesum.keys())],dtype=str).transpose()
        # sys.stderr.write("\nTotal effective %s pumping rates (Qi, Qf, Rate (s^-1) at 1 AU:\n" %(mol))
        # levels_lo = []
        # levels_up = []

        # for Qi,Qf in sorted(zip(QinitSum,QfinalSum)):
        #     level_lo = self.enlevels[self.enlevels['J'].str.match(Qi)]['Level'].values[0]
        #     level_up = self.enlevels[self.enlevels['J'].str.match(Qf)]['Level'].values[0]
        #     levels_lo.append(level_lo)
        #     levels_up.append(level_up)
        #     gratesumfinal[str(level_lo)+" "+str(level_up)]=gratesum[str(Qi)+" "+str(Qf)]
        # file_out = "g_"+str(self.mol)+"_1au.dat"
        # f = open(file_out,"w")
        # for Li,Lf in sorted(zip(levels_lo,levels_up)):
        #     f.write('%s %s %s \n' %(Li,Lf,gratesumfinal[str(Li)+" "+str(Lf)]))
        # f.close()




