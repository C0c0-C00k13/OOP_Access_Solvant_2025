"""Analyses of atom"""

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


def main():

    naccess_file, sasa_file = "Results/09-13-24/2c8r/2c8r.asa", "Results/06-21-2025/2c8r/SASA/output.asa"
    naccess_data = parse_atom_asa(naccess_file)
    sasa_data = parse_atom_asa(sasa_file)

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
    
    # print(all_keys)
    # print(x_labels)
    # print(y_values)
    # print(colors)
    # Display scatter plot
    display_atom_asa_diff(x_labels, y_values, colors)
