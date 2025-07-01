import numpy as np
from astroquery.linelists.cdms import CDMS
import astropy.units as u
import pandas as pd
import astropy.constants as const

def generate_Einstein_As(elo,lgint,freq,gup,Q):
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
    E_u = E_l + (freq.to(1/u.s) * const.h / const.k_B)

    intensity = 10**lgint * u.nm**2 * u.MHz

    gl_S_mu2 = 2.50251e4 * intensity.value * Q  / (freq.value * (np.exp(-E_l / (300*u.K)) - np.exp(-E_u / (300*u.K))))

    EA = 1.16395e-20 * (freq.value**3) * gl_S_mu2 / gup

    return EA, E_u.value

def getCDMS(mol, emax, min_frequency=0*u.GHz, max_frequency=5000*u.GHz, saveFile=True):
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
    ctbl = CDMS.query_lines(min_frequency=min_frequency, max_frequency=max_frequency, molecule=mol).to_pandas()

    #Drop NA columns to get rid of empty quantum numbers
    ctbl = ctbl.dropna(axis=1)

    #Drop transitions with Elo > Emax
    ctbl = ctbl[ctbl['ELO'] <= emax]

    #Drop 'name' and 'Lab' columns
    ctbl = ctbl.drop(['name','Lab'], axis=1)
    ctbl = ctbl.reset_index()

    #Get the molecule data
    molData = CDMS.get_species_table()
    mtbl = molData[molData['tag'] == int(mol.split()[0])]
    Q300 = 10**mtbl['lg(Q(300))'].item()

    #Calculate the Einstein A's and E_ups
    ctbl[['EA','EU']] = ctbl.apply(lambda x: pd.Series(generate_Einstein_As(x['ELO'],x['LGINT'],x['FREQ'],x['GUP'],Q300)), axis=1)

    #Work out how many quantum numbers we have. They come after QNFMT, which is column 9. 
    #Since we have upper and lower, divide by two to figure this out.
    #CDMS reserves at most 6 quanta for upper and lower levels
    cols = ctbl.columns
    qnfmt = np.where(cols == 'QNFMT')[0][0]
    nQuanta = int(len(cols[qnfmt+1:-2]))/2.
    if nQuanta == 1:
        ctbl['Qup'] = ctbl.apply(lambda x: str(int(x['Ju'])), axis=1)
        ctbl['Qlo'] = ctbl.apply(lambda x: str(int(x['Jl'])), axis=1)
    elif nQuanta == 2:
        ctbl['Qup'] = ctbl.apply(lambda x: str(int(x['Ju']))+'_'+str(int(x['Ku'])), axis=1)
        ctbl['Qlo'] = ctbl.apply(lambda x: str(int(x['Jl']))+'_'+str(int(x['Kl'])), axis=1)
    elif nQuanta == 3:
        ctbl['Qup'] = ctbl.apply(lambda x: str(int(x['Ju']))+'_'+str(int(x['Ku']))+'_'+str(int(x['vu'])), axis=1)
        ctbl['Qlo'] = ctbl.apply(lambda x: str(int(x['Jl']))+'_'+str(int(x['Kl']))+'_'+str(int(x['vl'])), axis=1)
    elif nQuanta == 4:
        ctbl['Qup'] = ctbl.apply(lambda x: str(int(x['Ju']))+'_'+str(int(x['Ku']))+'_'+str(int(x['vu']))+'_'+str(int(x['F1u'])), axis=1)
        ctbl['Qlo'] = ctbl.apply(lambda x: str(int(x['Jl']))+'_'+str(int(x['Kl']))+'_'+str(int(x['vl']))+'_'+str(int(x['F1l'])), axis=1)
    elif nQuanta == 5:
        ctbl['Qup'] = ctbl.apply(lambda x: str(int(x['Ju']))+'_'+str(int(x['Ku']))+'_'+str(int(x['vu']))+'_'+str(int(x['F1u']))+'_'+str(int(x['F2u'])), axis=1)
        ctbl['Qlo'] = ctbl.apply(lambda x: str(int(x['Jl']))+'_'+str(int(x['Kl']))+'_'+str(int(x['vl']))+'_'+str(int(x['F1l']))+'_'+str(int(x['F2l'])), axis=1)
    elif nQuanta == 6:
        ctbl['Qup'] = ctbl.apply(lambda x: str(int(x['Ju']))+'_'+str(int(x['Ku']))+'_'+str(int(x['vu']))+'_'+str(int(x['F1u']))+'_'+str(int(x['F2u']))+'_'+str(int(x['F3u'])), axis=1)
        ctbl['Qlo'] = ctbl.apply(lambda x: str(int(x['Jl']))+'_'+str(int(x['Kl']))+'_'+str(int(x['vl']))+'_'+str(int(x['F1l']))+'_'+str(int(x['F2l']))+'_'+str(int(x['F3l'])), axis=1)
    else:
        print('Problems defining quanta!')
        return

    #C34S and C33S are special cases where CDMS includes the v=0 and v=1 levels together
    if (mol.split()[1] == 'C34S') or (mol.split()[1] == 'C33S') or (mol.split()[1] == 'CS') or (mol.split()[1] == '13CS'):
        ctbl = ctbl[ctbl['Ku'] == 0]
        ctbl['Qup'] = ctbl.apply(lambda x: str(int(x['Ju'])), axis=1)
        ctbl['Qlo'] = ctbl.apply(lambda x: str(int(x['Jl'])), axis=1)

    #CN is a special case,too. it has v=0,1 and also needs to adjust the J,F QN's to half integer
    if mol.split()[1] == 'CN':
        ctbl = ctbl[ctbl['Ku'] == 0]
        ctbl['Qup'] = ctbl.apply(lambda x: str(int(x['Ju']))+'_'+str(int(x['vu'])-0.5)+'_'+str(int(x['F1u'])-0.5), axis=1)
        ctbl['Qlo'] = ctbl.apply(lambda x: str(int(x['Jl']))+'_'+str(int(x['vl'])-0.5)+'_'+str(int(x['F1l'])-0.5), axis=1)

    #Generate a list of levels and statistical weights
    enlevs = {}
    gs = {}
    ups = {}
    los = {}
    jlos = {}
    i=0
    for elo, gup, qup, qlo  in zip(ctbl['ELO'],ctbl['GUP'],ctbl['Qup'],ctbl['Qlo']):
        ups[i] = qup
        los[i] = qlo
        if mol.split()[1] == 'SO':
            print(qlo)
            jlos[qlo] = int(qlo.split('_')[1])
        else:
            jlos[qlo] = int(qlo.split('_')[0])
                        
        enlevs[qlo] = elo
        gs[qup] = gup
        i+=1

    #Work out number of energy levels
    nlev = len(enlevs.keys())
    #Start writing output if desired
    if saveFile:
        f = open(mol.split()[1]+'-out.dat', 'w')
        f.write("!MOLECULE\n")
        f.write(f"{mol.split()[1]:s}\n")
        f.write("!MOL. WEIGHT (AMU)\n")
        f.write(f"{ctbl['MOLWT'][0].item():d}\n")
        f.write("!NUMBER OF ENERGY LEVELS\n")
        f.write(f"{nlev:d}\n")
        f.write("!LEVEL + ENERGIES(cm^-1) + WEIGHT + J, F\n")

    sorted_enlevs = {key: value for key, value in sorted(enlevs.items(), key=lambda item: item[1])}
    enlevels = pd.DataFrame()
    levels = {}
    j=0
    for key in sorted_enlevs:
        j+=1

        #Lowest energy levels have undefined g, so set it according to J
        if key not in gs.keys():
            if mol.split()[1] == 'HCN':
                gs[key] = 3*(2*jlos[key]+1)
            else:
                gs[key] = (2*jlos[key] + 1)

        if j==1:
        #enlevels[j] = {}
            enlevels['Level'] = [j]
            enlevels['Energy'] = [enlevs[key]]
            enlevels['Weight'] = [gs[key]]
            enlevels['J'] = [str(key)]
        else:
            df = pd.DataFrame()
            df['Level'] = [j]
            df['Energy'] = [enlevs[key]]
            df['Weight'] = [gs[key]]
            df['J'] = [str(key)]

            enlevels = pd.concat([enlevels,df],ignore_index=True)

        levels[key] = j

        if saveFile:
            f.write(f"{j:5d}    {enlevs[key]:10.4f} {gs[key]:6.1f}   {key}\n")  

    if saveFile:
    #Count the number of valid transitions
        ntran = 0
        for l in np.arange(i):
            try:
                levels[ups[l]]
                ntran += 1
            except:
                pass
        f.write("!NUMBER OF RADIATIVE TRANSITIONS\n")
        f.write(f"{ntran:d}\n")
        f.write("!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)\n")

        tran=0
        for l in np.arange(i):
            try:
                levels[ups[l]]
                tran += 1
                f.write(f"{tran:5d} {levels[ups[l]]:5d} {levels[los[l]]:5d}   {ctbl['EA'][l]:10.4e}    {ctbl['FREQ'][l]/1000.:12.7f}   {ctbl['EU'][l]:7.2f}\n")
            except:
                pass

        f.write("!NUMBER OF COLL PARTNERS\n")
        f.write("1                        \n")
        f.write("!COLLISIONS BETWEEN      \n")                             
        f.write("1 - No collision data    \n")
        f.write("!NUMBER OF COLL TRANS    \n")
        f.write("0                        \n")
        f.write("!NUMBER OF COLL TEMPS    \n")
        f.write("0                        \n")
        f.write("!COLL TEMPS              \n")
        f.write("!NOTES:\n")
        f.write("!Generated using CDMS database and cdms2landa.py\n")

        f.close()

    return ctbl, enlevels