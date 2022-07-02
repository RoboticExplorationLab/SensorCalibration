# [src/MissionSim/CustomStructs/flags.jl]

""" FLAGS 

      Struct of flags that are used by the state machine to determine which state 
    to transition to.
"""
mutable struct FLAGS
    init_detumble::Bool 
    magnetometer_calibrated::Bool
    diodes_calibrated::Bool
    final_detumble::Bool 
    in_sun::Bool 

    function FLAGS(; in_sun = false, mag_cal = false, dio_cal = false, init_detumble = false, final_detumble = false)
        new(init_detumble, mag_cal, dio_cal, final_detumble, in_sun)
    end
end