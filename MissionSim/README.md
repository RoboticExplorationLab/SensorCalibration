# MISSION SIMULATOR
This repo contains code to run a complete "mission-style" simulation for the CubeSat project. Additional information can be found below or in the subfolders

## SYSTEM SETUP
During this simulation the satellite goes through four primary phases:
1. Magnetometer Calibration: The satellite calibrates the magnetometer on board (scale factor, biases, and non-orthogonality angles). Note that this method is attitude independent and allows us to correct future magnetometer readings. No controller is used during this phase.
2. Detumbling: The rotation of the satellite is slowed down by using a magnetorquer to generate torques opposite to the angular velocity. No state estimation is performed during this phase, and we hope that the generated magnetic fields do not affect the magnetometer calibration. 
3. Photodiode Calibration: The satellite calibrates the photodiodes on board (scale factor; azimuth and elevation angles of diodes) using an iterative non-linear least squares algorithm. No controller is used during this phase.
4. Magnetorquer Control: After calibration and detumbling, the satellite is ready to be controlled. A magnetorquer is used along with the Earth's magnetic field to orient the satellite as desired. A standard MEKF is used for attitude estimation. 

During each phase, the simulator iteratively updates the sensors, uses that data to update estimates, and then generates a command and uses the system dynamics to update the state until convergence. 