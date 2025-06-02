# accall.py
# VERSION: 2.1
# AIM:
# Input a Brookhaven entry file and output a PDB format file after filtering/cleaning,
# including Van der Waals radii, contained in an external file "vdw.radii".

# INPUT:
# PDB format file, van der Waals radii file

# OPTIONS:
# Inclusion of non-standard amino acids, het-atoms, waters, etc.
# Flagging of missing residues, chain-breaks, non-standard atom names, missing atoms.
# Nucleic acids, separable chains, polar/non-polar summing

# AUTHOR: S. Hubbard 3/92. EMBL.

import Naccess_Utils
import Accall_Params

# --- VARIABLES ---
i = j = k = 0
l = [0] * 256
ilen = flen = vlen = slen = nbackbone = 0

numats = [0] * maxr
resindex = [0] * maxs
backbone = [0] * maxs
polstats = [0] * maxs
rty = [0] * maxx
num_chains = num_res = nats = atype = nacids = 0
rtype = [0] * maxr
achain = [0] * maxs
rchain = [0] * maxx

spolar = [[0] * maxa for _ in range(maxr)]
vradii = [[0.0] * maxa for _ in range(maxr)]
xyz = [[0.0] * 3 for _ in range(maxs)]
rads = [0.0] * maxs
accs = [0.0] * maxs
occup = [0.0] * maxs
bfact = [0.0] * maxs
vdw = hyrad = 1.00
probe = zslice = 0.0

ressums = [[[0.0] * 2 for _ in range(7)] for _ in range(maxx)]
tsums = [0.0] * 7
csums = [0.0] * 7

# Characters
alt = firsta = ''
chnam = [''] * maxc
res = ''
aacids = [''] * maxr
rlab = ['RES', 'HEM', 'HOH']
atom = ''
anames = [[''] * maxa for _ in range(maxr)]
last = ''
resnam = [''] * maxx
label = [''] * maxs
fname = vname = sname = ''
card = ''
c = [''] * 256

# Logical flags
hetas = falt = hydro = fullo = conta = False
ok = start = wwaters = resok = aok = dorsa = oldr = False

# --- DEFAULTS ---
hetas = False
hydro = False
wwaters = False
fullo = False
dorsa = True
conta = False
oldr = False

# --- GET USER DIRECTIVES ---
while readstring(5, card, ilen) >= 0:
    n = parse(card, ilen, ' ', c, l)
    tolow(c[0], l[0])  # Assuming tolow modifies c[0] in place or returns a new value

    key = c[0][:4].lower()  # Fortran (1:4) → Python [:4]

    if key == 'pdbf':
        fname = c[1]
        flen = l[1]
    elif key == 'vdwf':
        vname = c[1]
        vlen = l[1]
    elif key == 'stdf':
        sname = c[1]
        slen = l[1]
    elif key == 'prob':
        probe = readfloat(c[1], l[1])
    elif key == 'zsli':
        zslice = readfloat(c[1], l[1])
    elif key == 'heta':
        hetas = True
    elif key == 'hydr':
        hydro = True
    elif key == 'wate':
        wwaters = True
    elif key == 'full':
        fullo = True
    elif key == 'oldr':
        oldr = True
    elif key == 'asao':
        dorsa = False
    elif key == 'cont':
        conta = True
    elif key == 'csid':
        nbackbone = 5

# --- OPEN FILES ---
i = fname.rfind('.')  # index of last '.' (Fortran: index(fname, '.') - 1)
k = 0
ok = False

# Find last '/' to isolate base filename
for j in range(i - 1, -1, -1):
    if fname[j] == '/' and not ok:
        k = j + 1
        ok = True

base = fname[k:i]  # filename without path and extension

if fopen(1, fname, flen, 'old') == 0:
    raise RuntimeError("ERROR: unable to open PDB file")

# Open output files
asa_filename = f"{base}.asa"
rsa_filename = f"{base}.rsa"
log_filename = f"{base}.log"

# Assuming Python-style file I/O or a wrapper around Fortran I/O
file2 = open(asa_filename, 'w')

if dorsa:
    file3 = open(rsa_filename, 'w')

file4 = open(log_filename, 'w')

# --- READ IN VDW RADII ---
vanin(
    vname,
    vlen,
    nacids,
    aacids,
    anames,
    numats,
    vradii,
    spolar,
    rtype
)

file4.write(' ACCALL - Accessibility calculations\n')
file4.write(f' MAX RESIDUES   {maxx:7d}\n')
file4.write(f' MAX ATOMS/RES  {maxa:7d}\n')
file4.write(f' PDB FILE INPUT {fname[:flen]}\n')
file4.write(f' PROBE SIZE     {probe:6.2f}\n')
file4.write(f' Z-SLICE WIDTH  {zslice:6.3f}\n')
file4.write(f' VDW RADII FILE {vname[:vlen]}\n')

file4.write(' INCL HETATOMS\n' if hetas else ' EXCL HETATOMS\n')
file4.write(' INCL HYDROGENS\n' if hydro else ' EXCL HYDROGENS\n')
file4.write(' INCL WATERS\n' if wwaters else ' EXCL WATERS\n')

file4.write(f' READDVW {nacids:3d} residues input\n')

# --- INITIALIZE ---
falt = False
start = True
last = ' ' * 10

nats = 0
num_res = 0
nchain = 1

# --- READ DATA & DECODE ---
while readstring(1, card, ilen) >= 0:
    atype = restype(card[:20])
    ip = -1

    if atype == 1 or (atype == 2 and hetas) or (atype == 3 and wwaters):
        # Ignore alternate positions other than blank or first encountered
        alt = card[16]
        if alt != ' ':
            if not falt:
                firsta = alt
                falt = True
            if alt != firsta:
                continue  # goto5

        # Ignore hydrogens & deuteriums unless flagged
        if card[13] in ('H', 'D', 'Q') or card[12] == 'H':
            if not hydro:
                continue  # goto5
            vdw = hyrad
            nats += 1
            continue  # goto6

        # Next atom
        nats += 1

        # First residue?
        if start:
            start = False
            chnam[0] = card[21]
            num_chains = 1

        # New residue?
        current_res = card[17:27]
        if last != current_res:
            last = current_res
            num_res += 1
            if num_res > maxx:
                file4.write(f"\n ERROR - Maximum number of residues exceeded {maxx}\n")
                file4.write(" Increase maxx in accall.pars and recompile\n")
                raise RuntimeError('SOLVA_ERROR: maxx exceeded')

            res = card[17:20]
            resnam[num_res - 1] = last

            i, resok = which3(res, aacids, nacids)
            ir = rtype[i] if i > 0 else 0
            if not resok:
                file4.write(f" WARNING - Unknown residue: {current_res}\n")

            rty[num_res - 1] = atype
            nchain = chain(num_chains, chnam, card[21])
            rchain[num_res - 1] = nchain

# --- Get atom type ---
atom = card[12:16]
backbone[nats] = what_atom(atom, ir, nbackbone)
achain[nats] = nchain

# --- Assign radius ---
ip = -1
if atom == ' OXT':
    vdw = 1.40
    ip = 1
else:
    if resok:
        vdw = 0.0
        i, aok, j = ratom(atom, anames, numats)
    else:
        aok = False

    if not resok or not aok:
        i, j, aok = gatom(atom, anames, nacids, numats)
        file4.write(f" WARNING - Unknown atom: {atom} in residue {card[17:27]}\n")
        if aok:
            file4.write(f" --> Using VDW {vradii[i][j]:.2f} for atom {atom} in {aacids[i]}\n")

    if not aok:
        vdw = vguess(atom)
        file4.write(f" --> GUESS VDW {vdw:.2f} for atom {atom} in {card[17:27]}\n")
    else:
        vdw = vradii[i][j]
        ip = spolar[i][j]

    if ip < 0:
        ip = polguess(atom)

# --- Store atom data ---
rads[nats] = vdw
polstats[nats] = ip
label[nats] = card[:30]
if fullo:
    occup[nats] = float(card[54:60])
    bfact[nats] = float(card[60:66])
resindex[nats] = num_res
xyz[nats] = list(map(float, [card[30:38], card[38:46], card[46:54]]))

# --- End of ATOM processing block ---
# 5 is the continue label in Fortran

# --- End of read loop ---
# End of: while readstring(1, card, ilen) >= 0

# --- Output ---
file1.close()
file4.write(' ADDED VDW RADII\n')
file4.write(f' CHAIN(S): {num_chains}, RESIDUES: {num_res}, ATOMS: {nats}\n')

# --- Calculate atomic accessibilities ---
accs = [0.0] * nats
solva(nats, xyz, rads, accs, probe, zslice)

if conta:
    for i in range(nats):
        accs[i] *= rads[i]**2 / (rads[i] + probe)**2

# --- Write ASA output ---
for i in range(nats):
    if fullo:
        file2.write(f"{label[i]:30}{xyz[i][0]:8.3f}{xyz[i][1]:8.3f}{xyz[i][2]:8.3f}"
                    f"{occup[i]:6.2f}{bfact[i]:5.1f}{accs[i]:7.3f} {rads[i]:5.2f}\n")
    else:
        file2.write(f"{label[i]:30}{xyz[i][0]:8.3f}{xyz[i][1]:8.3f}{xyz[i][2]:8.3f}"
                    f"{accs[i]:8.3f} {rads[i]:5.2f}\n")

file4.write(' CALCULATED ATOMIC ACCESSIBILITIES\n')

# --- Residue-level accessibility summary ---
if dorsa:
    summer(
        sname,
        slen,
        nats,
        accs,
        backbone,
        polstats,
        rtype,
        resindex,
        resnam,
        ressums,
        tsums
    )

    file4.write(' SUMMED ACCESSIBILITIES OVER RESIDUES\n')
    file3.write(header_format)  # Assuming defined format strings

    if oldr:
        file3.write(oldr_header_1)
        file3.write(oldr_header_2)
    else:
        file3.write(newr_header_1)
        file3.write(newr_header_2)

    for i in range(resindex[nats - 1]):
        if oldr:
            file3.write(res_format_old.format(
                rlab[rty[i]], resnam[i],
                *[ressums[i][j][0] for j in range(7)],
                *[ressums[i][j][1] for j in range(7)]
            ))
        else:
            file3.write(res_format_new.format(
                rlab[rty[i]], resnam[i],
                ressums[i][0][0], ressums[i][0][1],
                *[ressums[i][j][0] for j in range(3, 7)],
                *[ressums[i][j][1] for j in range(3, 7)]
            ))

    file3.write(chain_summary_header)

    for i in range(num_chains):
        if chnam[i] == ' ':
            chnam[i] = '_'
        csums = [0.0] * 7
        for j in range(num_res):
            if rchain[j] == i + 1:
                for k in range(7):
                    csums[k] += ressums[j][k][0]
        if oldr:
            file3.write(chain_format_old.format(i + 1, chnam[i], *csums))
        else:
            file3.write(chain_format_new.format(i + 1, chnam[i], csums[0], *csums[3:7]))

    if oldr:
        file3.write(total_format_old.format(*tsums[:7]))
    else:
        file3.write(total_format_new.format(tsums[0], *tsums[3:7]))

# --- FORMAT Definitions in Python ---

# 102: Unknown residue type
def format_102(resname):
    return f" UNKNOWN residue type.............> {resname:>10}"

# 104: Non-standard atom warning
def format_104(atom, resname):
    return f" NON-STANDARD atom.|{atom:<4}| in residue> {resname:>10}"

# 106: Assumed VDW radius
def format_106(atom, resname, vdw, ref_res):
    return f" ASSUMED vdw of {atom:<4} in {resname:>10} = {vdw:5.2f} (same as {ref_res:<3})"

# 108: Guessed VDW radius
def format_108(atom, resname, vdw):
    return f" GUESSED vdw of {atom:<4} in {resname:>10} = {vdw:5.2f}"

# 110: Summary of chains/residues/atoms
def format_110(chains, residues, atoms):
    return f" CHAINS   {chains:5d}\n RESIDUES {residues:5d}\n ATOMS    {atoms:5d}"

# 120: Header for summed accessibility file
def format_120(structure_name):
    return (f"REM  File of summed (Sum) and % (per.) "
            f"accessibilities for {structure_name}")

# 125 and 126: Column headers for residue-level accessibility (old/new format)
format_125 = ("REM RES _ NUM      All-atoms   Non-P-side   Polar-Side   "
              "Total-Side   Main-Chain    Non-polar    All polar")
format_126 = ("REM RES _ NUM      All-atoms   Total-Side   Main-Chain    "
              "Non-polar    All polar")

# 130 and 131: Column headers (absolute and relative values)
format_130 = ("REM                ABS   REL    ABS   REL    ABS   REL    "
              "ABS   REL    ABS   REL    ABS   REL    ABS   REL")
format_131 = ("REM                ABS   REL    ABS   REL    ABS   REL    "
              "ABS   REL    ABS   REL")

# 150: Residue-level values (old format)
def format_150(res_code, resname, values):
    return f"{res_code:<3} {resname:<10} " + " ".join(
        f"{v:7.2f}{p:6.1f}" for v, p in zip(values[::2], values[1::2]))

# 151: Residue-level values (new format)
def format_151(res_code, resname, values):
    return f"{res_code:<3} {resname:<10} " + " ".join(
        f"{v:7.2f}{p:6.1f}" for v, p in zip(values[::2], values[1::2]))

# 154: End of per-chain residue summaries
format_154 = "END  Absolute sums over single chains surface"

# 155: Chain-level sum (old format)
def format_155(chain_num, chain_name, csums):
    return f"CHAIN {chain_num:2d} {chain_name:1s}   " + "     ".join(
        f"{v:8.1f}" for v in csums)

# 156: Chain-level sum (new format)
def format_156(chain_num, chain_name, csums):
    return f"CHAIN {chain_num:2d} {chain_name:1s}   " + "     ".join(
        f"{v:8.1f}" for v in csums)

# 160: End of all-chain sums (old format)
def format_160(tsums):
    return "END  Absolute sums over all chains \nTOTAL        " + "     ".join(
        f"{v:8.1f}" for v in tsums)

# 161: End of all-chain sums (new format)
def format_161(tsums):
    return "END  Absolute sums over all chains \nTOTAL        " + "     ".join(
        f"{v:8.1f}" for v in tsums)
