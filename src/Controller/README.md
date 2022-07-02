# CONTROLLER
There are two primary controllers used: one to slow down the rotation of the satellite (detumbler) and one to orient the satellite. Both rely on a magnetorquer. 

All controller structs can be called with the `generate_command(ctrl)` command.

## Detumbler 
In order to slow the rotation of the satellite, the magnetorquer is used to generate a torque in opposition to the rotation of the satellite. There are a few different variations of this controller (B cross is the only one implemented; at some point, we may add B dot and a bang-bang style).

(*Magnetic Detumbling of a Rigid Spacecraft,* Avanzini 2012)

## Attitude Control 
TBD still