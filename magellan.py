# This code is a part of Maternal Genealogy Lineage Analyser - MaGelLAn-v2.0demo.
# MaGelLAn is an open source software and it is free for non-commercial use, as long as it is properly referenced.
# Authors are Ino Curik, Dalibor Hrsak and Strahil Ristov

import argparse
import sys
from mag_gui import runGUI
from mag_cl import MagellanCL

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('--mode', dest='mode', default= "gui", help='Magellan mode, GUI or CL')
    
    args = parser.parse_args()

    mode = args.mode
    if (mode != "gui") and (mode != "cl"):
        sys.stderr.write("Allowed values for --mode are gui or cl\n")
        sys.exit()

    return mode


def main():

    mode = parse_arguments()

    if mode == "gui":
        runGUI()
    elif mode == "cl":
        magellan_cl = MagellanCL()
        magellan_cl()

if __name__ == "__main__":
    main()
