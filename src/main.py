"""
Minimal working example of interval propagation through contact.

Run code with argument "simulate" to observe trajectory of a disk
with nominal parameters together with the lower and upper bound
trajectories from intervals in parameters.

Author: Nima Fazeli
lab: MCube
Adviser: Alberto Rodriguez
Department: Mechanical Engineering
Institution: MIT
Date: June 2018
"""

from classes import disk
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Minimal example for interval propagation through contact.')
    parser.add_argument('--simulate', help='Simulates a disk with nominal and bounded intervals.')
    parser.add_argument('--animation', help='Currently under construction, animates the simulations.')
    args = parser.parse_args()

    if args.simulate is None or args.simulate == 'simulate':
        # initialize the disk class and simulate it
        dsk = disk.Disk()
        dsk.simulate()

    if args.animation == 'animation':
        # run simulation and plot results
        dsk = disk.Disk()
        dsk.simulate()
        # work in progress, animating the results
        # dsk.animate_loop()

