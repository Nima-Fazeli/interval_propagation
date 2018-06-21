# Contact Model Interval Propagation
Code to propagate intervals of parameters and states through contact models.

## For a minimal example of a disk:

Run the following code:

'''python
cd src/
python main.py --simulate
'''

You'll get as output the trajectories of a disk with intervals defined over the parameters being dropped. The disk properties are:

'''python
R_range = [0.45, 0.55]
friction_range = [0.25, 0.35]
restitution_range = [0.45, 0.55]
mass_range = [0.95, 1.05]
'''

The initial conditions are:

'''python
q0 = [0.0, 0.7, 0.0]
v0 = [0.1, 0.0, 0.0]
'''

The trajectories are labeled as nominal, lower, and upper. The nominal is the middle value for the intervals of the parameters of the disk, and the lower and upper are the simulations of the two bounds of the trajectories of the disk, i.e. the two most distant possible trajectories.
