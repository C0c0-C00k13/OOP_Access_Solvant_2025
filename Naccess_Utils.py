import numpy as np
from math import acos, atan2, sqrt, pi

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
