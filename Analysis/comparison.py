""""""

import atomAnalysis
import residueAnalysis
import utils
import sys
import os

# Add parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(parent_dir)

# Import SASA.py
import SASA


if __name__ == "__main__":

    atomAnalysis.main()

    residueAnalysis.main()

    durations = utils.repeat_and_time(SASA.main, n=5)
    print(durations)



    # sasa_header = ["REM",  "CHAIN", "ResID", "ResName", "TotalASA", "RSA(%)", "MainASA", "MainRSA(%)", "SideASA", "SideRSA(%)", "PolarASA", "PolarRSA", "ApolarASA", "ApolarRSA"]
    # naccesse_header = ["REM", "RES", "CHAIN", "NUM", "All-atoms-ASA", "All-atoms-RSA", "Total-Side-ASA", "Total-Side-RSA", "Main-Chain-ASA", "Main-Chain-RSA", "Non-polar-ASA", "Non-polar-RSA", "All-polar-ASA", "All-polar-RSA"]

