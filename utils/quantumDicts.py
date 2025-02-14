def spin_to_int(level,mol):
    """ Converts spin species (e.g., A+,A-) to integers (+1,-1).
        This makes it easy to quantify a change between the two states
        via multiplication by -1. This is also useful for converting between
        LAMDA and HITRAN syntax

    """
    import re
    if mol == 'aCH3OH':
        m = re.match(' (.)(.) (.)(.)  A(.)     ',level)
        if m.group(5) == '-':
            quanta = m.group(1).rjust(2)+m.group(2)+m.group(3).rjust(2)+m.group(4)+'  -1     '
        else:
            quanta = m.group(1).rjust(2)+m.group(2)+m.group(3).rjust(2)+m.group(4)+'   1     '
    elif mol == 'eCH3OH':
        m = re.match(' (.)(.) (.)(.)  E(.)     ',level)
        if m.group(5) == '2':
            quanta = m.group(1).rjust(2)+m.group(2)+m.group(3).rjust(2)+m.group(4)+'  -1     '
        else:
            quanta = m.group(1).rjust(2)+m.group(2)+m.group(3).rjust(2)+m.group(4)+'   1     '
    return quanta

def generate_upper_quanta(level,mol):
    """
    Generate upper quanta for linear molecules in HITRAN
    """
    branch = level.split()[0]
    if 'e' in level.split()[1]:
        Jlo = level.split()[1].strip('e')
    elif 'f' in level.split()[1]:
        Jlo = level.split()[1].strip('f')
    elif 'F' in level.split()[1]:
        Jlo = level.split()[1].strip('F')
    elif 'E' in level.split()[1]:
        Jlo = level.split()[1].strip('E')
    else:
        Jlo = level.split()[1]
    Jlo = int(Jlo)
    if branch == 'P':
        Jup = str(abs(Jlo - 1))
    elif branch == 'R':
        Jup = str(Jlo + 1)
    elif branch == 'Q':
        Jup = str(Jlo)
    return Jup

def cleanup_lower_quanta(level,mol):
    """
    Reformat lower quanta for linear molecules in HITRAN
    """
    if 'e' in level.split()[1]:
        Jlo = level.split()[1].strip('e')
    elif 'f' in level.split()[1]:
        Jlo = level.split()[1].strip('f')
    elif 'F' in level.split()[1]:
        Jlo = level.split()[1].strip('F')
    elif 'E' in level.split()[1]:
        Jlo = level.split()[1].strip('E')
    else:
        Jlo = level.split()[1]
    Jlo = str(Jlo)
    return Jlo

def generate_statistical_weights(level,mol):
    """
    Generate statistical weights for HNC
    """
    if (mol == 'HNC'):
        J = int(level)
        g = 6*((2*J)+1)
    return g

def generate_ammonia_qns(level):
    import re
    s=re.match('(\d{1,2})\s{1,2}(\d{1,2}) ([as]) ([A][12].)([A][12].)',level)
    if s.group(4) == s.group(5):
        ns = '{:02d}_{:02d}_00'.format(int(s.group(1)),int(s.group(2)))
    else:
        ns = '{:02d}_{:02d}_01'.format(int(s.group(1)),int(s.group(2)))

    return ns