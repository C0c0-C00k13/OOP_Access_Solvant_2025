""""""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import io
import utils
from collections import defaultdict
import atomAnalysis
import residueAnalysis

# ----- INPUT DATA -----
# Replace these with actual file reads (e.g., open("file.txt").read())
naccess_data = """
RES GLY A   1    64.17  80.1  31.94  98.8  32.23  67.5  31.94  85.1  32.22  75.7
RES ILE A   2     0.64   0.4   0.06   0.0   0.58   1.6   0.06   0.0   0.58   1.6
RES VAL A   3    35.90  23.7  35.77  31.3   0.13   0.3  35.77  31.0   0.13   0.4
RES GLU A   4    85.07  49.4  81.01  60.1   4.06  10.8  51.28  85.1  33.79  30.2
RES THR B  27    81.09  58.2  69.85  68.7  11.24  29.9  62.22  82.2  18.87  29.7
RES PRO B  28   107.01  78.6  56.43  47.1  50.58 311.6  83.39  68.9  23.62 155.5
"""

tool_data = """
Chain ResID  ResName TotalASA   RSA(%)   PolarASA   PolarRSA   ApolarASA   ApolarRSA 
A     1      GLY     57.16      67.25    30.91      36.36      26.25       30.89     
A     2      ILE     2.38       1.31     2.38       1.31       0.00        0.00      
A     3      VAL     30.19      18.87    0.00       0.00       30.19       18.87     
A     4      GLU     81.31      44.43    41.93      22.91      39.38       21.52     
B     27     THR     83.53      57.21    24.46      16.75      59.07       40.46     
B     28     PRO     102.92     71.97    26.79      18.73      76.13       53.24     
"""


########################
## PARSE ATOM RESULTS
#######################

def display_difference_per_residue(tool_df):
    # ----- PARSE NACCESS OUTPUT -----
    naccess_pattern = r"RES (\w{3}) (\w) +(\d+) +([\d.]+) +([\d.]+)"
    matches = re.findall(naccess_pattern, naccess_data)
    naccess_df = pd.DataFrame(matches, columns=["ResName", "Chain", "ResID", "ASA", "RSA"])
    naccess_df["ResID"] = naccess_df["ResID"].astype(int)
    naccess_df["ASA"] = naccess_df["ASA"].astype(float)
    naccess_df["RSA"] = naccess_df["RSA"].astype(float)

    # ----- MERGE -----
    merged = pd.merge(tool_df, naccess_df, on=["Chain", "ResID", "ResName"], suffixes=('_tool', '_naccess'))

    # ----- CALCULATE DIFFERENCES -----
    merged["ASA_diff"] = merged["TotalASA"] - merged["ASA"]
    merged["RSA_diff"] = merged["RSA(%)"] - merged["RSA"]
    merged["ASA_abs_diff"] = merged["ASA_diff"].abs()
    merged["RSA_abs_diff"] = merged["RSA_diff"].abs()

    # ----- PLOTTING -----
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # ASA Plot
    sns.scatterplot(data=merged, x="ASA_abs_diff", y="ASA_diff", hue="ResName", ax=axes[0])
    axes[0].set_title("ASA Difference vs Absolute Difference")
    axes[0].set_xlabel("Absolute ASA Difference")
    axes[0].set_ylabel("ASA Difference")

    # RSA Plot
    sns.scatterplot(data=merged, x="RSA_abs_diff", y="RSA_diff", hue="ResName", ax=axes[1])
    axes[1].set_title("RSA Difference vs Absolute Difference")
    axes[1].set_xlabel("Absolute RSA Difference")
    axes[1].set_ylabel("RSA Difference")

    plt.tight_layout()
    plt.show()

# --------------------------------------------------




if __name__ == "__main__":
    # ----- PARSE TOOL OUTPUT -----
    tool_df = pd.read_csv(io.StringIO(tool_data), delim_whitespace=True)
    display_difference_per_residue(tool_df)
    
    atomAnalysis.main()
    
    residueAnalysis.main()

    

    sasa_header = ["REM",  "CHAIN", "ResID", "ResName", "TotalASA", "RSA(%)", "MainASA", "MainRSA(%)", "SideASA", "SideRSA(%)", "PolarASA", "PolarRSA", "ApolarASA", "ApolarRSA"]
    naccesse_header = ["REM", "RES", "CHAIN", "NUM", "All-atoms-ASA", "All-atoms-RSA", "Total-Side-ASA", "Total-Side-RSA", "Main-Chain-ASA", "Main-Chain-RSA", "Non-polar-ASA", "Non-polar-RSA", "All-polar-ASA", "All-polar-RSA"]
    
    # import time

    # def test_function():
    #     time.sleep(2)  # Simulate a task taking 2 seconds

    # durations = utils.repeat_and_time(test_function, n=3)

