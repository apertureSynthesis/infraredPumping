def omega(theta):
    import numpy as np
    """Solid angle subtended by the Sun
    Parameters
    ----------
    theta : float
        Angular diameter of the Sun
    Returns
    -------
    omega : float
        The resulting solid angle for the Sun
        https://en.wikipedia.org/wiki/Solid_angle
    """
    omega = 2*np.pi*(1.-np.cos(theta/2.))
    return omega


def glu(sigma, gu, gl, distance=1):
    import numpy as np
    import astropy.constants as const
    """Pumping rate from lower to upper level glu
    divided by the Einstein coefficient Aul
    (equation A7 from Crovisier & Encrenaz 1983)
    Parameters
    ----------
    sigma : float
        wavenumber in units of cm-1
    gu : float
        upper level degeneracy
    gl : float
        lower level degeneracy
    Returns
    -------
    glu : float
        This is actually glu/Aul
    """
    # solid angle subtended by the Sun
    om = omega(9.35e-3/distance)
    # black-body Sun temperature in K
    Tbb = 5778.
    # relative population of the lower level
    pl = 1.
    glu = (om/(4.*np.pi))*(gu/gl)*pl/(np.exp(const.h.value*const.c.value * sigma*1e2/const.k_B.value/Tbb)-1.)

    return glu

def generate_statistical_weights(level):
    """
    Generate statistical weights for HNC
    """
    J = int(level)
    g = 6*((2*J)+1)
    return g

def generate_linear_quanta(lower,upper,mol):

    """
    HITRAN doesn't include upper local quanta for linear molecules.
    Reformat the lower quanta into a more digestible format and calculate the upper quanta
    """
    import re
    
    #Dictionary of branch notation for linear molecules
    branch_dict = {'N': -3, 'O': -2, 'P': -1, 'Q': 0, 'R': 1, 'S': 2, 'T': 3}

    
    if mol in ['HCN', 'CO', 'CS', 'HNC', 'C34S', 'C33S', 'OCS', 'HC3N']:
        """
        CDMS format is simply J. HITRAN does not populate the upper quantum numbers at all, but lists the lower quantum numbers with P, Q, R notation.
        """
        lo = re.match('\s*([NOPSQRT])\s*(\d*)([a-z]*)(.*)',lower)
        loBranch = lo.group(1)
        Jlo = int(lo.group(2))
        Jup = int(Jlo + branch_dict[loBranch])

        lowerQuanta = f"{int(Jlo)}"
        upperQuanta = f"{int(Jup)}"

    elif mol in ['SO']:
        """
        CDMS format is N J. HITRAN does not populate the upper quantum numbers at all, but lists the lower quantum numbers with P, Q, R notation.
        """      
        lo = re.match('\s*([NOPQRST])\s*(\d*)([NOPQRST])\s*(\d*)(.*)(.*)',lower)  

        nBranch = lo.group(1)
        Nlo = int(lo.group(2))
        jBranch = lo.group(3)
        Jlo = int(lo.group(4))
        Nup = int(Nlo + branch_dict[nBranch])
        Jup = int(Jlo + branch_dict[jBranch])

        lowerQuanta = str(str(Nlo)+'_'+str(Jlo))
        upperQuanta = str(str(Nup)+'_'+str(Jup))

    else:
        lowerQuanta = f"{lower}"
        upperQuanta = f"{upper}"

    return lowerQuanta, upperQuanta

def generate_ammonia_qns(level):
    import re
    s=re.match('(\d{1,2})\s{1,2}(\d{1,2}) ([as]) ([A][12].)([A][12].)',level)
    if s.group(4) == s.group(5):
        ns = '{:02d}_{:02d}_00'.format(int(s.group(1)),int(s.group(2)))
    else:
        ns = '{:02d}_{:02d}_01'.format(int(s.group(1)),int(s.group(2)))

    return ns

def select_species(df,mol):
    """
    Split out a spin species from the HITRAN query if applicable.
    """
    if (mol == 'aCH3OH'):
        df = df[df['local_lower_quanta'].str.contains('A')]
    
    elif (mol == 'eCH3OH'):
        df = df[df['local_lower_quanta'].str.contains('E')]

    elif (mol == 'oH2O') or (mol == 'pH2O'):
        #For H2O, ortho (Ka-Kc=odd) vs para (Ka-Kc=even). Quantum numbers have format J_Ka_Kc
        df['is-ortho'] = df.apply(lambda x: True if (int(x['local_lower_quanta'].split()[1]) - int(x['local_lower_quanta'].split()[2])) %2 != 0 else False, axis=1)

        if (mol == 'oH2O'):
            df = df[df['is-ortho']]
        elif (mol == 'pH2O'):
            df = df[~df['is-ortho']]

    elif (mol == 'oH2CO') or (mol == 'pH2CO'):
        #For H2CO, ortho (Ka = odd) vs para (Ka = even). Quantum numbers have format J_Ka_Kc
        df['is-ortho'] = df.apply(lambda x: True if (int(x['local_lower_quanta'].split()[1] %2 != 0)) else False, axis=1)

        if (mol == 'oH2CO'):
            df = df[df['is-ortho']]
        elif (mol == 'pH2CO'):
            df = df[~df['is-ortho']]

    elif (mol == 'oNH3') or (mol == 'pNH3'):
        #For NH3, ortho (K = 3n) vs para (K = 3n+1, 3n+2). Quantum numbers have format J_K Sym.
        df['is-ortho'] = df.apply(lambda x: True if (int(x['local_lower_quanta'].split()[1] %3 == 0)) else False, axis=1)

        if (mol == 'oNH3'):
            df = df[df['is-ortho']]
        elif (mol == 'pNH3'):
            df = df[~df['is-ortho']]
    
    else:
        df = df

    return df

def hitran_to_cdms(lower,upper,mol):
    """
    Map the HITRAN local quanta to the CDMS local quanta format
    CDMS format is J, K, v, F1, F2, F3
    """
    import re
    
    #Dictionary of branch notation for linear molecules
    branch_dict = {'O': -2, 'P': -1, 'Q': 0, 'R': 1, 'S': 2, 'T': 3}

    
    if mol in ['HCN', 'CO', 'CS', 'HNC', 'C34S', 'C33S', 'OCS', 'HC3N']:
        """
        CDMS format is simply J. HITRAN does not populate the upper quantum numbers at all, but lists the lower quantum numbers with P, Q, R notation.
        We've already generated these into CDMS format in a previous step
        """

        lowerQuanta = str(lower)
        upperQuanta = str(upper)

    elif mol in ['SO']:
        """
        CDMS format is N J. HITRAN does not populate the upper quantum numbers at all, but lists the lower quantum numbers with P, Q, R notation.
        We've already generated these in CDMS format in a previous step.
        """      

        lowerQuanta = lower
        upperQuanta = upper

    elif mol in ['H2O','HDO','H2CO','SO2','H2S','HCOOH']:
        """
        CDMS format is J Ka Kc, HITRAN format is J Ka Kc
        """
        for level in [lower,upper]:
            lo = re.match('\s*(\d{1,3})\s*(\d{1,3})\s*(\d{1,3})(.{5})(.)',lower)
            up = re.match('\s*(\d{1,3})\s*(\d{1,3})\s*(\d{1,3})(.{5})(.)',upper)

            lowerQuanta = lo.group(1)+'_'+lo.group(2)+'_'+lo.group(3)
            upperQuanta = up.group(1)+'_'+up.group(2)+'_'+up.group(3)

    
    return lowerQuanta, upperQuanta


def hitran_to_lamda(lower,upper,mol):
    """
    Map the HITRAN local quanta to the LAMDA local quanta format
    """
    import re

    levels = []
    

    if mol in ['aCH3OH','eCH3OH']:
        """
        LAMDA format is J_(+/-)K, where the +/- is determined based on the spin species. HITRAN format is J K Spin.
        """
        for level in [lower,upper]:
            print(level)
            m = re.match('\s*(\d{1,2})\s*(\d{1,2})\s*([AE])([+-12])\s*',level)
            if m.group(4) == '-':
                #quanta = m.group(1).rjust(2)+m.group(2)+m.group(3).rjust(2)+m.group(4)+'  -1     '
                quanta = m.group(1).rjust(2)+'_-'+m.group(2)
            elif m.group(4) == '+':
                quanta = m.group(1).rjust(2)+'_'+m.group(2)

            elif m.group(4) == '2':
                quanta = m.group(1).rjust(2)+'_-'+m.group(2)
            elif m.group(4) == '1':
                quanta = m.group(1).rjust(2)+'_'+m.group(2)
            levels.append(quanta)

    elif mol in ['oH2O','pH2O','oH2CO','pH2CO']:
        """
        LAMDA format is J_Ka_Kc, HITRAN format is J Ka Kc
        """
        for level in [lower,upper]:
            m = re.match('\s*(\d{1,2})\s*(\d{1,2})\s*(\d{1,2})\s*',level)
            quanta = m.group(1)+'_'+m.group(2)+'_'+m.group(3)
            levels.append(quanta)

    elif mol in ['oNH3','pNH3']:
        """
        LAMDA format is J_K_S, where S is 0 if A or E state does not change and 1 if it does. HITRAN format is J K S.
        """
        for level in [lower,upper]:
            m=re.match('(\d{1,2})\s{1,2}(\d{1,2}) ([as]) ([AE][12].)([AE][12].)',level)
            if m.group(4) == m.group(5):
                quanta = '{:02d}_{:02d}_00'.format(int(m.group(1)),int(m.group(2)))
            else:
                quanta = '{:02d}_{:02d}_01'.format(int(m.group(1)),int(m.group(2)))
            levels.append(quanta)

    elif mol in ['HCN', 'CO', 'CS', 'HNC']:
        """
        LAMDA format is simply J. HITRAN does not populate the upper quantum numbers at all, but lists the lower quantum numbers with P, Q, R notation.
        """
        lowerQuanta = str(re.sub('([EFef])',' ',lower))
        cleanQuanta = int(re.sub('([PRQef]|1$)', ' ',lower,))
        levels.append(f"{cleanQuanta:02d}")

        branch = lowerQuanta.split()[0]
        Jlo = lowerQuanta.split()[1]
        Jlo = int(Jlo)
        if branch == 'P':
            upperQuanta = abs(Jlo - 1)
        elif branch == 'R':
            upperQuanta = abs(Jlo + 1)
        elif branch == 'Q':
            upperQuanta = abs(Jlo)
        levels.append(f"{upperQuanta:02d}")
    
    return levels[0], levels[1]

    


# def hitran_to_lamda(hdf,ldf,mol):
#     """
#     Map all local HITRAN quantum numbers to LAMDA format
#     """
#     #Generate upper quanta for HCN, CO, and CS since HITRAN doesn't list these
#     if (mol == 'HCN') or (mol == 'CS') or (mol == 'CO'):
#         hdf.loc[:,'Jf']=hdf['local_lower_quanta'].apply(generate_upper_quanta,args=(mol,))
#         hdf.loc[:,'Ji']=hdf['local_lower_quanta'].apply(cleanup_lower_quanta,args=(mol,))
        
#         ldf['J'] = ldf['J'].str.strip('"').map(int).map(str)

#     #Generate upper quanta and statistical weights for HNC since GEISA doesn't list these
#     if mol == 'HNC':
#         hdf.loc[:,'Jf']=hdf['local_lower_quanta'].apply(generate_upper_quanta,args=(mol,))
#         hdf.loc[:,'Ji']=hdf['local_lower_quanta'].apply(cleanup_lower_quanta,args=(mol,))
#         hdf.loc[:,'gp']=hdf['Jf'].apply(generate_statistical_weights,args=(mol,))
#         hdf.loc[:,'gpp']=hdf['Ji'].apply(generate_statistical_weights,args=(mol,))   

#         ldf['J'] = ldf['J'].str.strip('"').map(int).map(str)

#     #Split out ortho (K = 3n) from para (K = 3n+1, 3n+2) for NH3 from HITRAN

#     if mol == 'oNH3':
#         hdf['Jf'] = hdf['local_upper_quanta'].str.split().str[0]
#         hdf['Jf'] = hdf['Jf'].astype(int)
#         hdf['K0'] = hdf['local_lower_quanta'].str.split().str[1]
#         hdf['K0'] = hdf['K0'].astype(int)
#         hdf['J0'] = hdf['local_lower_quanta'].str.split().str[0]
#         hdf['J0'] = hdf['J0'].astype(int)
#         hdf = hdf[hdf['K0'] %3 == 0]
#         hdf = hdf[hdf['Jf'] <= 7]
#         hdf = hdf[hdf['J0'] <= 7]
#     elif mol == 'pNH3':
#         hdf['Jf'] = hdf['local_upper_quanta'].str.split().str[0]
#         hdf['Jf'] = hdf['Jf'].astype(int)
#         hdf['K0'] = hdf['local_lower_quanta'].str.split().str[1]
#         hdf['K0'] = hdf['K0'].astype(int)
#         hdf['J0'] = hdf['local_lower_quanta'].str.split().str[0]
#         hdf['J0'] = hdf['J0'].astype(int)
#         hdf = hdf[hdf['K0'] %3 != 0]
#         hdf = hdf[hdf['Jf'] <= 5]
#         hdf = hdf[hdf['J0'] <= 5]   
    
#     #Split out ortho (Ka = odd) vs para (Ka = even) for H2CO from HITRAN
#     if mol == 'oH2CO':
#         hdf['J0'],hdf['Ka0'],hdf['Kc0'] = hdf['local_lower_quanta'].str.split().str
#         hdf['Ka0'] = hdf['Ka0'].astype(int)
#         hdf = hdf[hdf['Ka0'] % 2 != 0]
#     elif mol == 'pH2CO':
#         hdf['J0'],hdf['Ka0'],hdf['Kc0'] = hdf['local_lower_quanta'].str.split().str
#         hdf['Ka0'] = hdf['Ka0'].astype(int)
#         hdf = hdf[hdf['Ka0'] % 2 == 0]

#     #Split out ortho (Ka-Kc = odd) vs para (Ka-Kc = even) for H2O from HITRAN
#     if mol == 'oH2O':
#         hdf['J0'],hdf['Ka0'],hdf['Kc0'] = hdf['local_lower_quanta'].str.split().str
#         hdf['Ka0'] = hdf['Ka0'].astype(int)
#         hdf['Kc0'] = hdf['Kc0'].astype(int)
#         hdf = hdf[(hdf['Ka0']-hdf['Kc0']) %2 != 0]
#     elif mol == 'pH2O':
#         hdf['J0'],hdf['Ka0'],hdf['Kc0'] = hdf['local_lower_quanta'].str.split().str
#         hdf['Ka0'] = hdf['Ka0'].astype(int)
#         hdf['Kc0'] = hdf['Kc0'].astype(int)
#         hdf = hdf[(hdf['Ka0']-hdf['Kc0']) %2 == 0]

#     #Convert CH3OH spin states from letters to integers from HITRAN: A+ --> +1, A- --> -1, E1 --> +1, E2 --> -1.
#     #Then find Qinit, Qup, and Qfinal for CH3OH. I've called the quantum numbers "Q" since for nonlinear molecules there's more than just J
#     if (mol == 'aCH3OH') or (mol =='eCH3OH'):
#         #Split out either A or E methanol
#         if (mol == 'aCH3OH'):
#             species='A'
#         else:
#             species='E'
#         hdf = hdf['lower_lower_quanta'].str.contains(species)
#         #Pull out the lower and upper quanta for each transition
#         hdf.loc[:,'local_upper_quanta_int']=hdf['local_upper_quanta'].apply(spin_to_int,args=(mol,))
#         hdf.loc[:,'local_lower_quanta_int']=hdf['local_lower_quanta'].apply(spin_to_int,args=(mol,))