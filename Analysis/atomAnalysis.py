"""Analyses of atom"""

import utils
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Optional for hover tooltips
try:
    import mplcursors
except ImportError:
    mplcursors = None


def parse_atom_asa(file_path):
    """
    Parses a file and extracts atom number, atom name, and ASA.
    Returns a dictionary with keys as "atom_num-atom_name" and values as ASA (float).
    """
    asa_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                parts = line.split()
                atom_num = parts[1]
                atom_name = parts[2]
                asa = float(parts[-2])
                key = f"{atom_num}-{atom_name}"
                asa_dict[key] = asa
    return asa_dict

def display_atom_asa_diff(x_labels, y_values, colors):
    """
    Displays ASA values per atom from different sources, with labels next to each point.
    """
    plt.figure(figsize=(12, 6))
    scatter = plt.scatter(range(len(x_labels)), y_values, c=colors, s=30)

    # Add text labels next to each point
    for i, label in enumerate(x_labels):
        plt.text(i + 0.05, y_values[i], label, fontsize=6, rotation=45, color='black')

    plt.xlabel('Atom Index')
    plt.ylabel('ASA Value')
    plt.title('ASA per Atom from NACCESS and SASA.py')
    plt.xticks([])  # Hide tick labels on x-axis since labels are next to points
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    if mplcursors:
        cursor = mplcursors.cursor(scatter, hover=True)
        @cursor.connect("add")
        def on_add(sel):
            sel.annotation.set_text(f"{x_labels[sel.index]}\nASA: {y_values[sel.index]:.3f}")
    else:
        print("Install mplcursors for hover interaction: pip install mplcursors")

    # Add legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='NACCESS', markerfacecolor='blue', markersize=8),
        Line2D([0], [0], marker='o', color='w', label='SASA.py', markerfacecolor='red', markersize=8)
    ]
    plt.legend(handles=legend_elements, title="Data Source", loc="upper left")

    plt.show()

def display_naccess_vs_sasa(naccess_data, sasa_data):
    """
    Plots NACCESS ASA (x-axis) vs SASA.py ASA (y-axis) per atom.
    Labels each point and distinguishes missing data with color.
    """
    x_vals = []
    y_vals = []
    labels = []
    colors = []

    all_keys = set(naccess_data) | set(sasa_data)
    for key in sorted(all_keys, key=lambda x: int(x.split('-')[0])):  # sort by atom number
        nacc = naccess_data.get(key)
        sasa = sasa_data.get(key)

        # Only plot if both values are present
        if nacc is not None and sasa is not None:
            x_vals.append(nacc)
            y_vals.append(sasa)
            labels.append(key)
            colors.append('purple')
        elif nacc is not None:
            x_vals.append(nacc)
            y_vals.append(0)
            labels.append(key)
            colors.append('blue')
        elif sasa is not None:
            x_vals.append(0)
            y_vals.append(sasa)
            labels.append(key)
            colors.append('red')

    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(x_vals, y_vals, c=colors, s=30)

    for i, label in enumerate(labels):
        plt.text(x_vals[i] + 0.3, y_vals[i], label, fontsize=6, color='black', rotation=30)

    plt.xlabel("ASA (NACCESS)")
    plt.ylabel("ASA (SASA.py)")
    plt.title("NACCESS vs SASA.py ASA Values per Atom")
    plt.grid(True, linestyle='--', alpha=0.6)

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Both Present', markerfacecolor='purple', markersize=6),
        Line2D([0], [0], marker='o', color='w', label='NACCESS Only', markerfacecolor='blue', markersize=6),
        Line2D([0], [0], marker='o', color='w', label='SASA.py Only', markerfacecolor='red', markersize=6),
    ]
    plt.legend(handles=legend_elements, title="Data Availability", fontsize=8, title_fontsize=9)

    plt.tight_layout()
    plt.show()

def main():

    PDB_IDS = {'1':'1bj5', '2':'1c26', '3' : '2c8r', '4' : '2oe4', '5' : '6pwf'}
    response = utils.timed_input("Choose a PDB ID:\n\
1':'1bj5'\n'2':'1c26'\n'3' : '2c8r'\n'4' : '2oe4'\n'5' : '6pwf'\n",10)
    pdb_id = PDB_IDS.get(response,'2c8r')
    naccess_file = f"Results/NACCESS/{pdb_id}/{pdb_id}.asa"
    sasa_file = f"Results/SASA/06-24-2025/{pdb_id}/output.asa"
    naccess_data = parse_atom_asa(naccess_file)
    sasa_data = parse_atom_asa(sasa_file)

    response = utils.timed_input("Select \n1 : naccess & sasa on the same axis\n\
2 : naccess vs sasa\nDefault : naccess vs sasa\n", 7)
    if response == '1':
        # Merge keys and prepare lists
        all_keys = set(naccess_data) | set(sasa_data)
        x_labels = []
        y_values = []
        colors = []

        for key in all_keys:
            if key in naccess_data:
                x_labels.append(key)
                y_values.append(naccess_data[key])
                colors.append('blue')
            if key in sasa_data:
                x_labels.append(key)
                y_values.append(sasa_data[key])
                colors.append('red')

        # Display scatter plot
        display_atom_asa_diff(x_labels, y_values, colors)
    else:
        display_naccess_vs_sasa(naccess_data, sasa_data)

if __name__ == "__main__":
    main()
