# SIMULATOR
This folder contains a simulator used to propogate state dynamics and generate the measurements that the various filters use. Because this file essentially plays the role of the environment, it will not be used in the actual flight and is merely for generating data to work with. 

## Simulator Module 
This module exports 'rk4', the primary function for updating state, and 'generate_measurements', which takes in the state and generates 
both idealized and noisy measurements.

Example function calls:   
``` next_state = rk4(inertia_matrix, state, control, time, time_step)```   
``` truth, sensors, eclipse_factor, noise = generate_measurements(satellite_struct, albedo_struct, state, time, time_step)```   

