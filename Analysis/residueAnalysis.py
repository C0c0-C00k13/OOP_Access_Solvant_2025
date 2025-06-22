""""""
import matplotlib.pyplot as plt
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
                        total_asa = float(parts[4])  # column 4 (index 4)
                        asa[(chain, resnum)] = total_asa
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
                        total_asa = float(parts[4])  # column 3 (index 4)
                        asa[(chain, resnum)] = total_asa
                    except ValueError:
                        continue
    return asa

def plot_total_asa(naccess_asa, sasapy_asa):
    """
    Plot total ASA from NACCESS vs SASA.py for residues present in both.
    """
    x_vals = []
    y_vals = []
    labels = []

    # Use residues present in both
    common_keys = set(naccess_asa.keys()) & set(sasapy_asa.keys())
    for key in sorted(common_keys):
        x = naccess_asa[key]
        y = sasapy_asa[key]
        x_vals.append(x)
        y_vals.append(y)
        labels.append(f"{key[0]}{key[1]}")

    plt.figure(figsize=(8, 8))
    plt.scatter(x_vals, y_vals, color='blue', alpha=0.7)
    plt.xlabel("Total ASA (NACCESS)")
    plt.ylabel("Total ASA (SASA.py)")
    plt.title("Total ASA Comparison per Residue")
    plt.grid(True)

    # Optional: annotate points with residue label
    for i, label in enumerate(labels):
        plt.annotate(label, (x_vals[i], y_vals[i]), textcoords="offset points", xytext=(3,3), ha='left', fontsize=8)

    plt.tight_layout()
    plt.show()

def main():
    file_naccess, file_py = "Results/09-13-24/2c8r/2c8r.rsa", "Results/06-22-2025/2c8r/SASA/output.rsa"
    
    naccess_asa = parse_naccess_total_asa(file_naccess)
    sasapy_asa = parse_sasapy_total_asa(file_py)
    print(naccess_asa)
    print(sasapy_asa)

    plot_total_asa(naccess_asa, sasapy_asa)