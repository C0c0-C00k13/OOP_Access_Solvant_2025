import threading
import time
import argparse


def parse_args():
    """Parses the arguments provided in the command line.

    Returns
    ---
    parser.parse_args() (Namespace) : Argument inputs.
    """

    parser = argparse.ArgumentParser(description="Compute ASA/RSA from PDB.")

    parser.add_argument("pdb_file", type=str, help="PDB file to process")
    parser.add_argument("-i", "--hetero", choices=["y", "n"], default="n",
                        help="Include HETATM records (y/n), default: n")
    parser.add_argument("-n", "--points", type=int, default=92,
                        help=f"Number of points representing the sphere, default: {92}")
    parser.add_argument("-o", "--output", type=str, default="output",
                        help="Name of the output file, default : output")
    parser.add_argument("-p", "--probe", type=float, default=1.4,
                        help="Probe radius (default: 1.4)")
    parser.add_argument("-r", "--radii", type=str, default=None,
                        help="Custom radii file (default: None)")

    return parser.parse_args()

def timed_input(prompt, timeout=5):
    """"""
    user_input = [None]

    def get_input():
        user_input[0] = input(prompt)

    thread = threading.Thread(target=get_input)
    thread.daemon = True
    thread.start()
    thread.join(timeout)

    if thread.is_alive():
        print("\nTime expired!")
        return None
    else:
        return user_input[0]

def repeat_and_time(func, n=1, *args, **kwargs):
    durations = []

    for i in range(1, n + 1):
        print(f"\nRun {i}/{n}...")
        start_time = time.time()

        # Run the target function
        func(*args, **kwargs)

        end_time = time.time()
        elapsed = end_time - start_time
        durations.append(elapsed)

        mins, secs = divmod(elapsed, 60)
        print(f"Duration: {int(mins)} min {secs:.2f} sec")

    return durations