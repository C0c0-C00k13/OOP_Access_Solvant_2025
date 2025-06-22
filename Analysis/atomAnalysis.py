"""Analyses of atom"""

import matplotlib.pyplot as plt
import utils
from collections import defaultdict

########################
## PARSE ATOM RESULTS
#######################

def parse_files(file1:str,file2:str):
    with open(file1) as f1, open(file2) as f2:
        file1_lines = f1.readlines()
        file2_lines = f2.readlines()

    timelimit = 10
    response = utils.timed_input(f"Enter something in {timelimit} seconds:\n\
avg : average ASA of each atom_name\n\
default : for ASA of each atom\t", timelimit)
    if response.lower() == "avg" :
        differences = mean_column10_by_atom(file1_lines, file2_lines)
        x_labels = differences.keys()
        diff = differences.values()
    else:
        differences = compare_files_by_atom_and_column(file1_lines, file2_lines)
        x_labels = [atom_diff[0] for atom_diff in differences]
        diff = [atom_diff[1] for atom_diff in differences]
    display_atom_asa_diff(x_labels, diff)


def compare_files_by_atom_and_column(file1_lines, file2_lines):
    result = []
    
    for line1, line2 in zip(file1_lines, file2_lines):
        # Skip lines that aren't ATOM records
        if not (line1.startswith("ATOM") and line2.startswith("ATOM")):
            continue
        
        columns1 = line1.split()
        columns2 = line2.split()

        try:
            atom = columns1[1]  
            y1, y2 = float(columns1[9]), float(columns2[9])
            diff = abs(y1 - y2)
            result.append((atom, diff))
        except (IndexError, ValueError):
            continue  # Skip malformed lines

    return result


def mean_column10_by_atom(file1_lines, file2_lines):
    atom_values = defaultdict(list)

    for line1, line2 in zip(file1_lines, file2_lines):
        if not (line1.startswith("ATOM") and line2.startswith("ATOM")):
            continue

        columns1 = line1.split()
        columns2 = line2.split()

        try:
            atom = columns1[2]  # 2nd column (atom name)
            y1 = float(columns1[9])
            y2 = float(columns2[9])
            atom_values[atom].extend([y1, y2])
        except (IndexError, ValueError):
            continue  # skip malformed lines

    # Calculate mean for each atom type
    mean_by_atom = {
        atom: sum(values) / len(values)
        for atom, values in atom_values.items()
    }

    return mean_by_atom

########################
## GRAPH
#######################

def display_atom_asa_diff(x_labels, differences):
    # Create bar plot
    plt.figure(figsize=(10, 6))
    # bars = plt.bar(x_labels, differences, color='skyblue')
    scatters = plt.scatter(x_labels, differences, color='skyblue')
    plt.xlabel('Atom (with index)')
    plt.ylabel('Absolute Difference in ASA')
    plt.title('Differences ASA Between NACCESS and SASA.py')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add hover interactivity
    try:
        import mplcursors
        # cursor = mplcursors.cursor(bars, hover=True)
        cursor = mplcursors.cursor(scatters, hover=True)
        @cursor.connect("add")
        def on_add(sel):
            sel.annotation.set_text(f"{x_labels[sel.index]}, {differences[sel.index]:.3f}")
    except ImportError:
        print("The 'mplcursors' module is not installed. You can install it with:")
        print("pip install mplcursors")
        mplcursors = None  # Optional: Set to None to check later

        
    # Show the plot
    plt.show()

def main():
    file_naccess, file_py = "Results/09-13-24/2c8r/2c8r.asa", "Results/06-21-2025/2c8r/SASA/output.asa"
    parse_files(file_naccess, file_py)