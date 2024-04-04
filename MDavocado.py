#!/usr/bin/env python3
from ramachandran import *
import sys

def guacamole(topology,trajectory):
    """
    The main wrapper function that initialises the RamachandranPlots class
    from the MDAnalysis Universe and makes all MDavocado plots for the given trajectory.
    """
    u = Universe(topology,trajectory)
    rp = RamachandranPlots(u)
    rp.run()

if __name__ == "__main__":
    #TODO make this much nicer and friendlier with argparse library
    topology, trajectory = sys.argv[1], sys.argv[2]
    guacamole(topology,trajectory)
