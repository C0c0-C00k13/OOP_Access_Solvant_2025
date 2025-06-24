""""""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import utils


def parse_naccess_total_asa(file_naccess):
    """
    Parse NACCESS file lines to extract total ASA per residue.
    Returns dict keyed by (chain, resnum) with total ASA float value.
    """
    asa = {}
    with open(file_naccess, 'r') as f1:
        for line in f1:
            if line.startswith("RES"):
                parts = line.split()
                if len(parts) >= 5:
                    chain = parts[2]
                    resnum = parts[3]
                    try:
                        total_asa = float(parts[4]) # Total atoms
                        side_asa = parts[6] # Total-Side
                        main_asa = parts[8] # Total-Main
                        apolar_asa = parts[10] # Non-polar
                        polar_asa = parts[12] # All polar
                        asa[(chain, resnum)] = {"total_asa" : total_asa,
                                                "side_asa" : side_asa,
                                                "main_asa" : main_asa,
                                                "apolar_asa" : apolar_asa,
                                                "polar_asa" : polar_asa}
                    except ValueError:
                        continue
    return asa

def parse_sasapy_total_asa(file_py):
    """
    Parse SASA.py file lines to extract total ASA per residue.
    Returns dict keyed by (chain, resnum) with total ASA float value.
    """
    asa = {}
    with open(file_py, 'r') as f2:
        for line in f2:
            if line.startswith("RES"):
                parts = line.strip().split()
                # print(line)
                # print(parts)
                if len(parts) >= 6:
                    chain = parts[1]
                    resnum = parts[2]
                    try:
                        total_asa = float(parts[4]) # Total atoms
                        main_asa = parts[6] # Main
                        side_asa = parts[8] # Side
                        polar_asa = parts[10] # Polar
                        apolar_asa = parts[12] # Apolar
                        asa[(chain, resnum)] = {"total_asa" : total_asa,
                                                "main_asa" : main_asa,
                                                "side_asa" : side_asa,
                                                "polar_asa" : polar_asa,
                                                "apolar_asa" : apolar_asa}
                    except ValueError:
                        continue
    return asa


import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D

def plot_asa_comparison_by_type(naccess_asa, sasapy_asa):
    """
    Show both NACCESS and SASA.py ASA values on the same scatterplot.
    Points are colored based on data source (NACCESS or SASA.py).
    """
    valid_types = {'1':'total_asa', '2':'main_asa','3': 'side_asa','4': 'polar_asa','5' :'apolar_asa'}
    response = utils.timed_input("Choose a type of ASA comparison:\n\
1 : 'total_asa'\n2 : 'main_asa'\n3 : 'side_asa'\n\
4 : 'polar_asa'\n5 : 'apolar_asa'\nDefault : total_asa\n", timeout=20)
    
    asa_type = valid_types.get(response, 'total_asa')

    x_vals = []
    y_vals = []
    colors = []
    labels = []

    # Source-specific color
    source_color = {
        'NACCESS': 'blue',
        'SASA.py': 'red'
    }

    # Use common residues
    common_keys = sorted(set(naccess_asa.keys()) & set(sasapy_asa.keys()))
    residue_indices = {key: idx for idx, key in enumerate(common_keys)}

    for key in common_keys:
        idx = residue_indices[key]
        label = f"{key[0]}{key[1]}"

        try:
            naccess_val = float(naccess_asa[key][asa_type])
            x_vals.append(idx)
            y_vals.append(naccess_val)
            colors.append(source_color['NACCESS'])
            labels.append(label)
        except (KeyError, ValueError):
            continue

        try:
            sasapy_val = float(sasapy_asa[key][asa_type])
            x_vals.append(idx)
            y_vals.append(sasapy_val)
            colors.append(source_color['SASA.py'])
            labels.append(label)
        except (KeyError, ValueError):
            continue

    plt.figure(figsize=(10, 6))
    plt.scatter(x_vals, y_vals, c=colors, alpha=0.7)

    # Annotate every residue once (at the NACCESS point)
    for i in range(0, len(x_vals), 2):  # step by 2 to only annotate once per residue
        plt.annotate(labels[i], (x_vals[i], y_vals[i]), textcoords="offset points", xytext=(3, 3), fontsize=8)

    plt.xticks(range(len(common_keys)), [f"{k[0]}{k[1]}" for k in common_keys], rotation=45)
    plt.ylabel(f"{asa_type.replace('_', ' ').title()} Value")
    plt.title(f"{asa_type.replace('_', ' ').title()} – NACCESS vs SASA.py")
    plt.grid(True)

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='NACCESS', markerfacecolor='blue', markersize=8),
        Line2D([0], [0], marker='o', color='w', label='SASA.py', markerfacecolor='red', markersize=8)
    ]
    plt.legend(handles=legend_elements, title="Data Source", loc="upper right")

    plt.tight_layout()
    plt.show()



def main():
    PDB_IDS = {'1':'1bj5', '2':'1c26', '3' : '2c8r', '4' : '2oe4', '5' : '6pwf'}
    response = utils.timed_input("Choose a PDB ID:\n\
1':'1bj5'\n'2':'1c26'\n'3' : '2c8r'\n'4' : '2oe4'\n'5' : '6pwf'\n",10)
    pdb_id = PDB_IDS.get(response,'2c8r')
    file_naccess = f"Results/NACCESS/{pdb_id}/{pdb_id}.rsa"
    file_py = f"Results/SASA/06-24-2025/{pdb_id}/output.rsa"

    naccess_asa = parse_naccess_total_asa(file_naccess)
    sasapy_asa = parse_sasapy_total_asa(file_py)

    plot_asa_comparison_by_type(naccess_asa, sasapy_asa)

if __name__ == "__main__":
    main()
