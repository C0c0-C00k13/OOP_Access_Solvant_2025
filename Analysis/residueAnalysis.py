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


def plot_asa_comparison_by_type(naccess_asa, sasapy_asa):
    """
    Plot ASA comparison for a specific ASA type between NACCESS and SASA.py.
    Color points based on ASA type.
    
    asa_type: 'total_asa', 'main_asa', 'side_asa', 'polar_asa', or 'apolar_asa'
    """
    valid_types = {'1':'total_asa', '2':'main_asa','3': 'side_asa','4': 'polar_asa','5' :'apolar_asa'}
    response = utils.timed_input("Choose a type ofcomparison asa\n\
1 : 'total_asa'\n\
2 : 'main_asa'\n\
3 : 'side_asa'\n\
4 : 'polar_asa'\n\
5 : 'apolar_asa'\n", timeout=20)
    
    asa_type = valid_types.get(response, 'total_asa')
    x_vals = []
    y_vals = []
    colors = []
    labels = []

    # Define a color map for each ASA type
    asa_color_map = {
        'total_asa': 'blue',
        'main_asa': 'green',
        'side_asa': 'orange',
        'polar_asa': 'purple',
        'apolar_asa': 'red'
    }

    common_keys = set(naccess_asa.keys()) & set(sasapy_asa.keys())
    for key in sorted(common_keys):
        try:
            x = float(naccess_asa[key][asa_type])
            y = float(sasapy_asa[key][asa_type])
            x_vals.append(x)
            y_vals.append(y)
            colors.append(asa_color_map[asa_type])
            labels.append(f"{key[0]}{key[1]}")
        except (KeyError, ValueError):
            continue  # Skip if data is missing or not a float

    plt.figure(figsize=(8, 8))
    scatter = plt.scatter(x_vals, y_vals, color=colors, alpha=0.7, label=asa_type)
    plt.xlabel(f"{asa_type.replace('_', ' ').title()} (NACCESS)")
    plt.ylabel(f"{asa_type.replace('_', ' ').title()} (SASA.py)")
    plt.title(f"{asa_type.replace('_', ' ').title()} Comparison per Residue")
    plt.grid(True)

    # Annotate each point
    for i, label in enumerate(labels):
        plt.annotate(label, (x_vals[i], y_vals[i]), textcoords="offset points", xytext=(3, 3), ha='left', fontsize=8)

    # Add legend for ASA type
    plt.legend(title="ASA Type")
    plt.tight_layout()
    plt.show()


def main():
    file_naccess, file_py = "Results/09-13-24/2c8r/2c8r.rsa", "Results/06-22-2025/2c8r/SASA/output.rsa"
    
    naccess_asa = parse_naccess_total_asa(file_naccess)
    sasapy_asa = parse_sasapy_total_asa(file_py)
    # print(naccess_asa)
    # print(sasapy_asa)

    plot_asa_comparison_by_type(naccess_asa, sasapy_asa)