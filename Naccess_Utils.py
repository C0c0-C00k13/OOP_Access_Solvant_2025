import numpy as np
from math import acos, atan2, sqrt, pi

def vguess(atom: str) -> float:
    """
    Guess the van der Waals radius for a given atom name.

    Parameters:
        atom (str): A 4-character string representing the atom name.

    Returns:
        float: Guessed van der Waals radius.
    """
    atom = atom.ljust(4)  # Ensure it's at least 4 characters
    vdw = 1.80  # default guess

    if atom[1] == 'C':
        vdw = 1.80
    if atom[1] == 'N':
        vdw = 1.60
    if atom[1] == 'S':
        vdw = 1.85
    if atom[1] == 'O':
        vdw = 1.40
    if atom[1] == 'P':
        vdw = 1.90
    if atom[:2] == 'CA':
        vdw = 2.07
    if atom[:2] == 'FE':
        vdw = 1.47
    if atom[:2] == 'CU':
        vdw = 1.78
    if atom[:2] == 'ZN':
        vdw = 1.39
    if atom[:2] == 'MG':
        vdw = 1.73

    return vdw


def chain(c, names):
    """
    Return integer number of chain assigned from the single-letter ID.
    If not present, assign a new index and add to names.

    Parameters:
    - c (str): A single character representing the chain ID.
    - names (list): A list of existing chain names (single-character strings).

    Returns:
    - chain_idx (int): The 1-based index of the chain in names.
    """
    for i, name in enumerate(names):
        if c == name:
            return i + 1  # 1-based indexing as in Fortran
    names.append(c)
    return len(names)


def gatom(atom, anames, nres, nats):
    """
    Find residue and atom index for a given atom name.

    Parameters:
    - atom (str): Atom name (4 characters).
    - anames (list of list of str): anames[residue][atom] - 2D list of atom names.
    - nres (int): Number of residues.
    - nats (list of int): Number of atoms in each residue.

    Returns:
    - ok (bool): True if found.
    - ir (int): Residue index (1-based).
    - ia (int): Atom index (1-based within residue).
    """
    ok = False
    ir = 0
    while ir < nres and not ok:
        ir += 1
        ok, ia = ratom(atom, anames, ir, nats)
    if not ok:
        ir = 0
        ia = 0
    return ok, ir, ia


def ratom(atom, anames, ires, nats):
    """
    Check if a standard atom name exists in a residue.

    Parameters:
    - atom (str): Atom name (4 characters).
    - anames (list of list of str): anames[residue][atom] - 2D list of atom names.
    - ires (int): Residue index (1-based).
    - nats (list of int): Number of atoms per residue.

    Returns:
    - ok (bool): True if atom found.
    - find (int): Atom index (1-based) if found, else 0.
    """
    ok = False
    find = 0
    if ires == 0 or ires > len(nats):
        return ok, find

    for i in range(nats[ires - 1]):
        if atom == anames[ires - 1][i]:
            find = i + 1  # 1-based index
            ok = True
            return ok, find

    return ok, find


def restype(card):
    """
    Determine residue type from a PDB record line.

    Parameters:
    - card (str): PDB line (at least 20 characters).

    Returns:
    - int: 
        1 if line starts with 'ATOM',
        2 if line starts with 'HETATM',
        3 if residue name at positions 18–20 is 'HOH',
        0 otherwise.
    """
    if card[0:4] == 'ATOM':
        return 1
    if card[0:6] == 'HETATM':
        return 2
    if card[17:20] == 'HOH':
        return 3
    return 0


def which3(res, acids, nacids):
    """
    Search for a 3-letter residue name in the list of known acids.

    Parameters:
    - res (str): 3-letter residue name to search.
    - acids (list of str): List of known 3-letter residue names.
    - nacids (int): Number of residues currently in acids list.

    Returns:
    - tuple: (ires, ok)
        ires (int): Index of found residue (1-based). If not found, returns nacids + 1.
        ok (bool): True if residue found, False otherwise.
    """
    for i in range(nacids):
        if res == acids[i]:
            return i + 1, True  # Fortran is 1-based
    return nacids + 1, False


def fopen(iochan, filename, flen, fstat):
    """
    Simulate opening a file with Fortran-style interface.

    Parameters:
    - iochan (int): Simulated file unit (not used here, placeholder).
    - filename (str): Filename to open.
    - flen (int): Length of the filename to use.
    - fstat (str): File status ("old", "new", "unknown", etc.).

    Returns:
    - int: 1 if file opened successfully, 0 if error occurred.
    """
    try:
        mode = {
            "old": "r",
            "new": "w",
            "unknown": "a+"
        }.get(fstat.lower(), "r")

        with open(filename[:flen], mode):
            pass
        return 1
    except Exception:
        return 0


def readstring(file_obj):
    """
    Read a line from a file and trim trailing spaces.

    Parameters:
    - file_obj: Open file object (Python file object, not an integer).

    Returns:
    - flen (int): Length of the trimmed line, or -1 on error or EOF.
    - card (str): The full line read from the file (up to 256 characters).
    """
    try:
        card = file_obj.readline()
        if not card:
            return -1, ""
        card = card[:256]
        flen = len(card.rstrip())
        return flen, card
    except Exception:
        return -1, ""


def parse(card, separator):
    """
    Parse a string into fields separated by a given separator.

    Parameters:
    - card (str): The input string.
    - separator (str): A single character used as the separator.

    Returns:
    - parse_count (int): Number of parsed elements.
    - chars (list of str): The parsed substrings.
    - clen (list of int): Lengths of each parsed substring.
    """
    chars = []
    clen = []

    length = len(card)
    i = 0
    search = False
    start = 0

    while i < length:
        if not search:
            if card[i] != separator:
                start = i
                search = True
        else:
            if card[i] == separator:
                chars.append(card[start:i])
                clen.append(i - start)
                search = False
        i += 1

    if search:
        chars.append(card[start:length])
        clen.append(length - start)

    parse_count = len(chars)
    return parse_count, chars, clen


def readint(card, clen):
    """
    Tries to parse an integer from a substring of `card` up to length `clen`.

    Returns:
        int: The parsed integer, or -999 on error.
    """
    try:
        value = float(card[:clen].strip())
        return int(value)
    except ValueError:
        return -999


def readfloat(card, clen):
    """
    Tries to parse a float from a substring of `card` up to length `clen`.

    Returns:
        float: The parsed float, or -999.9 on error.
    """
    try:
        return float(card[:clen].strip())
    except ValueError:
        return -999.9


def tolow(text, clen):
    """
    Converts the first `clen` characters of `text` to lowercase.

    Args:
        text (str): The input string.
        clen (int): Number of characters to convert.

    Returns:
        str: Modified string with lowercase letters up to `clen`.
    """
    return text[:clen].lower() + text[clen:]


def solva(nats, xyz, rads, probe, zslice, nac, ncube, nint):
    maxs = len(rads)

    accs = np.zeros(maxs)
    cube = np.zeros(maxs, dtype=int)
    rad = np.zeros(maxs)
    radsq = np.zeros(maxs)
    dx = np.zeros(nint)
    dy = np.zeros(nint)
    d = np.zeros(nint)
    dsq = np.zeros(nint)
    arci = np.zeros(nint)
    arcf = np.zeros(nint)
    tag = np.zeros(nint, dtype=int)
    inov = np.zeros(nint, dtype=int)

    natm = np.zeros((nac, ncube), dtype=int)
    itab = np.zeros(ncube, dtype=int)

    xmin, ymin, zmin = [9999.0] * 3
    xmax, ymax, zmax = [-9999.0] * 3

    # Setup constants
    pix2 = 2.0 * pi
    ict = nint

    # Preprocess: radii + probe, bounding box
    rmax = 0.0
    for i in range(nats):
        rad[i] = rads[i] + probe
        radsq[i] = rad[i] ** 2
        rmax = max(rmax, rad[i])
        xmin = min(xmin, xyz[i][0])
        ymin = min(ymin, xyz[i][1])
        zmin = min(zmin, xyz[i][2])
        xmax = max(xmax, xyz[i][0])
        ymax = max(ymax, xyz[i][1])
        zmax = max(zmax, xyz[i][2])

    rmax *= 2.0

    # Grid dimensions
    idim = int((xmax - xmin) / rmax + 1)
    if idim < 3: idim = 3
    jidim = int((ymax - ymin) / rmax + 1)
    if jidim < 3: jidim = 3
    jidim *= idim
    kjidim = int((zmax - zmin) / rmax + 1)
    if kjidim < 3: kjidim = 3
    kjidim *= jidim

    if kjidim > ncube:
        raise RuntimeError("SOLVA_ERROR: max cubes exceeded")

    # Fill cubes
    for l in range(nats):
        i = int((xyz[l][0] - xmin) / rmax + 1)
        j = int((xyz[l][1] - ymin) / rmax)
        k = int((xyz[l][2] - zmin) / rmax)
        kji = k * jidim + j * idim + i
        n = itab[kji] + 1
        if n > nac:
            raise RuntimeError("SOLVA_ERROR: max atoms per cube exceeded")
        itab[kji] = n
        natm[n - 1][kji] = l
        cube[l] = kji

    # Main loop over atoms
    nzp = int(1.0 / zslice + 0.5)

    for ir in range(nats):
        kji = cube[ir]
        io = 0
        area = 0.0
        xr, yr, zr = xyz[ir]
        rr = rad[ir]
        rrx2 = rr * 2.0
        rrsq = radsq[ir]

        # Find neighboring atoms in adjacent cubes
        for k in [-1, 0, 1]:
            for j in [-1, 0, 1]:
                for i in [-1, 0, 1]:
                    mkji = kji + k * jidim + j * idim + i
                    if 0 <= mkji < kjidim:
                        nm = itab[mkji]
                        if nm >= 1:
                            for m in range(nm):
                                in_atom = natm[m][mkji]
                                if in_atom != ir:
                                    if io >= ict:
                                        raise RuntimeError("SOLVA_ERROR: intrsctns > max")
                                    dx[io] = xr - xyz[in_atom][0]
                                    dy[io] = yr - xyz[in_atom][1]
                                    dsq[io] = dx[io] ** 2 + dy[io] ** 2
                                    d[io] = sqrt(dsq[io])
                                    inov[io] = in_atom
                                    io += 1

        if io >= 1:
            zres = rrx2 / nzp
            zgrid = xyz[ir][2] - rr - zres / 2
        else:
            accs[ir] = pix2 * rrx2 * rr
            continue

        for _ in range(nzp):
            zgrid += zres
            rsec2r = rrsq - (zgrid - zr) ** 2
            if rsec2r <= 0.0:
                continue
            rsecr = sqrt(rsec2r)
            karc = 0

            for j in range(io):
                in_atom = inov[j]
                rsec2n = radsq[in_atom] - (zgrid - xyz[in_atom][2]) ** 2
                if rsec2n <= 0.0:
                    continue
                rsecn = sqrt(rsec2n)
                if d[j] >= rsecr + rsecn:
                    continue
                b = rsecr - rsecn
                if d[j] <= abs(b):
                    continue

                trig_test = (dsq[j] + rsec2r - rsec2n) / (2.0 * d[j] * rsecr)
                trig_test = max(min(trig_test, 0.99999), -0.99999)
                alpha = acos(trig_test)
                beta = atan2(dy[j], dx[j]) + pi
                ti = beta - alpha
                tf = beta + alpha
                if ti < 0.0:
                    ti += pix2
                if tf > pix2:
                    tf -= pix2
                arci[karc] = ti
                if tf < ti:
                    arcf[karc] = pix2
                    karc += 1
                    arci[karc] = 0.0
                arcf[karc] = tf
                karc += 1

            if karc == 0:
                arcsum = pix2
            else:
                # Sort arcs
                indices = np.argsort(arci[:karc])
                arcsum = arci[indices[0]]
                t = arcf[indices[0]]
                for k in indices[1:]:
                    if t < arci[k]:
                        arcsum += arci[k] - t
                    t = max(t, arcf[k])
                arcsum += pix2 - t

            parea = arcsum * zres
            area += parea

        accs[ir] = area * rr

    print("SOLVA: PROGRAM ENDS CORRECTLY")
    return accs


def sortag(a):
    """
    Sorts the list `a` in ascending order and returns the original indices in `tag`.
    This replicates the behavior of the Fortran `sortag` subroutine.
    
    Parameters:
    a (list of float): The array to sort.

    Returns:
    tuple: (sorted_a, tag) where `sorted_a` is the sorted version of `a` and
           `tag` is a list of the original indices corresponding to the sorted elements.
    """
    n = len(a)
    tag = list(range(n))
    a_tagged = list(zip(a, tag))


    sorted_tagged = quicksort_with_tag(a_tagged)
    sorted_a = [x[0] for x in sorted_tagged]
    sorted_tag = [x[1] for x in sorted_tagged]
    
    return sorted_a, sorted_tag


def quicksort_with_tag(arr):
    if len(arr) <= 1:
        return arr
    pivot = arr[len(arr) // 2][0]
    left = [x for x in arr if x[0] < pivot]
    middle = [x for x in arr if x[0] == pivot]
    right = [x for x in arr if x[0] > pivot]
    return quicksort_with_tag(left) + middle + quicksort_with_tag(right)


def summer(
    sname,
    nats,
    accs,
    backbone,
    polstats,
    rtype,
    resindex,
    resnam,
    ressums,
    tsums,
    readstring,
    fopen,
    maxr,
    maxs,
    maxx
):
    """
    Summarize atomic accessibilities by residue.
    
    Parameters are assumed to be consistent with Fortran-style data
    (e.g., `accs`, `backbone`, `resindex`, etc. as lists of length `nats`).
    """
    stand = False
    standarea = [[0.0] * 7 for _ in range(maxr)]
    acids = [''] * maxr
    rindex = [0] * maxx

    # Try to open the standard accessibility file
    if fopen(1, sname, len(sname), 'old') != 0:
        print(f'REM  Relative accessibilites read from external file "{sname}"', file=open(3, 'a'))
        stand = True
        i = 0
        while True:
            line, ilen = readstring(1)
            if line is None or i >= maxr:
                break
            if line[0:4] == 'ATOM':
                acids[i] = line[12:15]
                standarea[i][0] = float(line[16:23])
                standarea[i][1] = float(line[29:36])
                standarea[i][2] = float(line[42:49])
                standarea[i][3] = float(line[55:62])
                standarea[i][4] = float(line[68:75])
                standarea[i][5] = float(line[81:88])
                standarea[i][6] = float(line[94:108])
                i += 1
        with open(4, 'a') as f:
            f.write(f' RELATIVE (STANDARD) ACCESSIBILITIES READFOR {i:3d} AMINO ACIDS\n')
    else:
        with open(4, 'a') as f:
            f.write(' NO STANDARD VALUES INPUT\n')

    nacids = i

    for i in range(resindex[nats - 1]):
        rindex[i] = 0
        if stand:
            res = resnam[i][:3]
            ok, ires = which3(res, acids, nacids)
            if ok:
                rindex[i] = ires

    # Initialize tsums and ressums if not already
    for i in range(len(tsums)):
        tsums[i] = 0.0

    for i in range(nats):
        ir = resindex[i]
        acc = accs[i]
        tsums[0] += acc
        ressums[ir][0][0] += acc

        if backbone[i] == 0:
            ressums[ir][4][0] += acc
            tsums[4] += acc
        else:
            ressums[ir][3][0] += acc
            tsums[3] += acc
            if polstats[i] == 0:
                ressums[ir][1][0] += acc
                tsums[1] += acc
            elif polstats[i] == 1:
                ressums[ir][2][0] += acc
                tsums[2] += acc

        if polstats[i] == 0:
            ressums[ir][5][0] += acc
            tsums[5] += acc
        else:
            ressums[ir][6][0] += acc
            tsums[6] += acc

    for i in range(resindex[nats - 1]):
        ires = rindex[i]
        if stand and ires != 0:
            for j in range(7):
                if standarea[ires][j] > 0.0:
                    ressums[i][j][1] = 100.0 * ressums[i][j][0] / standarea[ires][j]
                else:
                    ressums[i][j][1] = 0.0
        else:
            for j in range(7):
                ressums[i][j][1] = -99.9


def polguess(atom):
    """
    Guess polarity based on the second character of the atom name.

    Args:
        atom (str): Atom name, expected to be at least 2 characters long.

    Returns:
        int: 1 if the atom is polar (second char is 'O', 'N', or 'A'), else 0.
    """
    if len(atom) >= 2 and atom[1] in ('O', 'N', 'A'):
        return 1
    return 0


def what_atom(atom: str, ir: int, n: int) -> int:
    """
    Determine if the given atom matches standard main chain or nucleotide atoms.

    Args:
        atom (str): Atom name (4 characters expected, padding if needed).
        ir (int): Residue type indicator (1 for protein, 2 for nucleotide, others for general).
        n (int): Number of atoms in main chain (used for ir == 1 case).

    Returns:
        int: 0 if atom is recognized as standard, 1 otherwise.
    """
    mc = [' N  ', ' C  ', ' O  ', ' OXT', ' CA ']
    nc = [' P  ', ' O1P', ' O2P', ' O5*', ' C5*', ' C4*',
          ' O4*', ' C3*', ' O3*', ' C2*', ' C1*']

    if ir == 1:
        for i in range(n):
            if atom == mc[i]:
                return 0
    elif ir == 2:
        for i in range(11):
            if atom == nc[i]:
                return 0
    else:
        for i in range(4):
            if atom == mc[i]:
                return 0
        for i in range(11):
            if atom == nc[i]:
                return 0

    return 1


def vanin(
    vname: str,
    vlen: int,
    nacids: list,
    aacids: list,
    anames: list,
    numats: list,
    vradii: list,
    spolar: list,
    rtype: list,
    maxr: int,
    maxa: int
):
    try:
        with open(vname[:vlen], 'r') as f:
            nacids_val = 0
            for line in f:
                card = line.rstrip('\n')
                ilen = len(card.rstrip())
                n, c, l = parse(card, ilen, ' ')
                if c[0] == 'RESIDUE':
                    nacids_val += 1
                    if nacids_val > maxr:
                        raise RuntimeError("ERROR: increase maxr")
                    rtype_val = 1
                    if c[1][:4] == 'NUCL':
                        rtype_val = 2
                    elif c[1][:4] == 'HETA':
                        rtype_val = 3
                    rtype[nacids_val - 1] = rtype_val

                    aa3 = c[2][:3].replace('_', ' ')
                    aacids[nacids_val - 1] = aa3
                    numats[nacids_val - 1] = 0

                elif c[0] == 'ATOM':
                    idx = nacids_val - 1
                    numats[idx] += 1
                    if numats[idx] > maxa:
                        raise RuntimeError("ERROR: increase maxa")

                    atom_name = card[5:9].replace('_', ' ')
                    anames[idx][numats[idx] - 1] = atom_name

                    vrad = readfloat(card[10:14], 4)
                    vradii[idx][numats[idx] - 1] = vrad

                    if n >= 4:
                        pol = readint(card[15:16], 1)
                    else:
                        pol = -1
                    spolar[idx][numats[idx] - 1] = pol

            nacids[0] = nacids_val

    except FileNotFoundError:
        raise RuntimeError('ERROR: unable to open "vdw.radii"')

