# CellModelREU
# Authors: Rob Lorch, Alan Gan
# Advised by Dr. Barber, using his code/ideas as a resource
# Basic guide on how to use:

To run a simulation quickly, just use one of the scripts. 

To set up your own simulation, first call 'initializeNetwork', passing in the number of external nodes as well as the simulation type. Then, deform the cell with a call to 'deformCellForce' or 'deformCellDisplacement'. Finally, run the simulation with a call to an ODE solver (stepEuler, ode15s, etc.). For more detail, you can use one of the scripts as a template and/or review function documentation.

NOTE: Currently, the time integrator (using ode15s or stepEuler) is much quicker and more robust than the steady state solver (findSteadyState function).

