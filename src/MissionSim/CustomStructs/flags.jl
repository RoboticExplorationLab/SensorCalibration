# [src/MissionSim/CustomStructs/flags.jl]

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