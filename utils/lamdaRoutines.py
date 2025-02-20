def collider_ids():
    """
    List of possible collisional partner IDs in a LAMDA datafile
    """
    collider_ids = {'H2': 1,
                'PH2': 2,
                'OH2': 3,
                'E': 4,
                'H': 5,
                'HE': 6,
                'H+': 7}
    collider_ids.update({v: k for k, v in list(collider_ids.items())})
    return collider_ids

def parse_lamda_lines(data):
    """
    Extract a LAMDA datafile into a dictionary of tables. Credit Brian Svoboda

    (non-pythonic!  more like, fortranic)
    """

    meta_rad = {}
    meta_mol = {}
    meta_coll = {}
    levels = []
    radtrans = []
    collider = None
    ncolltrans = None
    for ii, line in enumerate(data):
        if line[0] == '!':
            continue
        if 'molecule' not in meta_mol:
            meta_mol['molecule'] = _cln(line)
            continue
        if 'molwt' not in meta_mol:
            meta_mol['molwt'] = float(_cln(line))
            continue
        if 'nenergylevels' not in meta_mol:
            meta_mol['nenergylevels'] = int(_cln(line))
            continue
        if len(levels) < meta_mol['nenergylevels']:
            lev, en, wt = _cln(line).split()[:3]
            jul = " ".join(_cln(line).split()[3:])
            levels.append([int(lev), float(en), int(float(wt)), jul])
            continue
        if 'radtrans' not in meta_rad:
            meta_rad['radtrans'] = int(_cln(line))
            continue
        if len(radtrans) < meta_rad['radtrans']:
            # Can have wavenumber at the end.  Ignore that.
            trans, up, low, aval, freq, eu = _cln(line).split()[:6]
            radtrans.append([int(trans), int(up), int(low), float(aval),
                             float(freq), float(eu)])
            continue
        if 'ncoll' not in meta_coll:
            meta_coll['ncoll'] = int(_cln(line))
            collrates = {}
            continue
        if collider is None:
            collider = int(line[0])
            collname = collider_ids[collider]
            collrates[collider] = []
            meta_coll[collname] = {'collider': collname,
                                   'collider_id': collider}
            continue
        if ncolltrans is None:
            ncolltrans = int(_cln(line))
            meta_coll[collname]['ntrans'] = ncolltrans
            continue
        if 'ntemp' not in meta_coll[collname]:
            meta_coll[collname]['ntemp'] = int(_cln(line))
            continue
        if 'temperatures' not in meta_coll[collname]:
            meta_coll[collname]['temperatures'] = [int(float(x)) for x in
                                                   _cln(line).split()]
            continue
        if len(collrates[collider]) < meta_coll[collname]['ntrans']:
            trans, up, low = [int(x) for x in _cln(line).split()[:3]]
            temperatures = [float(x) for x in _cln(line).split()[3:]]
            collrates[collider].append([trans, up, low] + temperatures)
        if len(collrates[collider]) == meta_coll[collname]['ntrans']:
            # meta_coll[collider_ids[collider]+'_collrates'] = collrates
            log.debug("{ii} Finished loading collider {0:d}: "
                      "{1}".format(collider, collider_ids[collider], ii=ii))
            collider = None
            ncolltrans = None
            if len(collrates) == meta_coll['ncoll']:
                # All done!
                break

    if len(levels[0]) == 4:
        mol_table_names = ['Level', 'Energy', 'Weight', 'J']
    elif len(levels[0]) == 5:
        mol_table_names = ['Level', 'Energy', 'Weight', 'J', 'F']
    else:
        raise ValueError("Unrecognized levels structure.")
    mol_table_columns = [table.Column(name=name, data=data)
                         for name, data in zip(mol_table_names,
                                               zip(*levels))]
    mol_table = table.Table(data=mol_table_columns, meta=meta_mol)

    rad_table_names = ['Transition', 'Upper', 'Lower', 'EinsteinA',
                       'Frequency', 'E_u(K)']
    rad_table_columns = [table.Column(name=name, data=data)
                         for name, data in zip(rad_table_names,
                                               zip(*radtrans))]
    rad_table = table.Table(data=rad_table_columns, meta=meta_rad)

    coll_tables = {collider_ids[collider]: None for collider in collrates}
    for collider in collrates:
        collname = collider_ids[collider]
        coll_table_names = (['Transition', 'Upper', 'Lower'] +
                            ['C_ij(T={0:d})'.format(tem) for tem in
                             meta_coll[collname]["temperatures"]])
        coll_table_columns = [table.Column(name=name, data=data)
                              for name, data in zip(coll_table_names,
                                                    zip(*collrates[collider]))]
        coll_table = table.Table(data=coll_table_columns,
                                 meta=meta_coll[collname])
        coll_tables[collname] = coll_table

    return coll_tables, rad_table, mol_table

def _cln(s):
    """
    Clean a string of comments, newlines
    """
    return s.split("!")[0].strip()

def write_lamda_datafile(filename, tables):
    """
    Write tuple of tables with LAMDA data into a datafile that follows the
    format adopted for the LAMDA database. Credit Brian Svoboda

    Parameters
    ----------
    filename : str
        Fully qualified path of the file to write.

    tables: tuple
        Tuple of Tables ({rateid: coll_table}, rad_table, mol_table)
    """
    import platform
    import sys

    collrates, radtransitions, enlevels = tables

    levels_hdr = ("""! MOLECULE
                  {0}
                  ! MOLECULAR WEIGHT
                  {1}
                  ! NUMBER OF ENERGY LEVELS
                  {2}
                  ! LEVEL + ENERGIES(cm^-1) + WEIGHT + J
                  """)
    levels_hdr = re.sub('^ +', '', levels_hdr, flags=re.MULTILINE)
    radtrans_hdr = ("""! NUMBER OF RADIATIVE TRANSITIONS
                    {0}
                    ! TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)
                    """)
    radtrans_hdr = re.sub('^ +', '', radtrans_hdr, flags=re.MULTILINE)
    coll_hdr = ("""! NUMBER OF COLL PARTNERS
                {0}
                """)
    coll_hdr = re.sub('^ +', '', coll_hdr, flags=re.MULTILINE)
    coll_part_hdr = ("""! COLLISION PARTNER
                     {0} {1}
                     ! NUMBER OF COLLISIONAL TRANSITIONS
                     {2}
                     ! NUMBER OF COLLISION TEMPERATURES
                     {3}
                     ! COLLISION TEMPERATURES
                     {4}
                     ! TRANS + UP + LOW + RATE COEFFS(cm^3 s^-1)
                     """)
    coll_part_hdr = re.sub('^ +', '', coll_part_hdr, flags=re.MULTILINE)

    if platform.system() == 'Windows':
        if sys.version_info[0] >= 3:
            stream = open(filename, 'w', newline='')
        else:
            stream = open(filename, 'wb')
    else:
        stream = open(filename, 'w')

    with stream as f:
        f.write(levels_hdr.format(enlevels.meta['molecule'],
                                  enlevels.meta['molwt'],
                                  enlevels.meta['nenergylevels']))
        enlevels.write(f, format='ascii.no_header')
        f.write(radtrans_hdr.format(radtransitions.meta['radtrans']))
        radtransitions.write(f, format='ascii.no_header')
        f.write(coll_hdr.format(len(collrates)))
        for k, v in collrates.items():
            temperatures = ' '.join([str(i) for i in v.meta['temperatures']])
            f.write(coll_part_hdr.format(v.meta['collider_id'],
                                         v.meta['collider'],
                                         v.meta['ntrans'],
                                         v.meta['ntemp'],
                                         temperatures))
            v.write(f, format='ascii.no_header')

def get_lamda_file(lamdamol):
    """
    Function to retrieve the latest line of molecular line files from the LAMDA datbase
    """
    from bs4 import BeautifulSoup
    import requests

    #Location of the LAMDA data files
    main_url = 'http://home.strw.leidenuniv.nl/~moldata/'

    #Get all the possible files
    soup = BeautifulSoup(requests.get(main_url).content, "html.parser")

    #Create a lookup dictionary from the links
    mol_dict = dict([(i,j) for i,j in zip((link["href"].split('.')[0] for link in soup.select('a[href*=".dat"]')), (main_url+link["href"] for link in soup.select('a[href*=".dat"]')))])

    #Download the requested file and save it locally
    response = requests.get(mol_dict[lamdamol])

    with open(lamdamol+".dat", "w") as f:
        f.write(response.text)


def query_lamda(lamdamol):
    import os
    """
    Query a file from the LAMDA database, parse it, and put it into a dataframe
    """

    #Check if the LAMDA file has already been retrieved
    pwd = os.getcwd()
    if not os.path.exists(pwd+'/'+lamdamol+'.dat'):
        get_lamda_file(lamdamol)
    
    f = open(pwd+'/'+lamdamol+'.dat')
    datafile = [s.strip() for s in f]
    collrates, radtransitions, enlevels = parse_lamda_lines(datafile)

    return collrates, radtransitions, enlevels

