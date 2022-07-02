# [src/MissionSim/CustomStructs/noise.jl]

"""
    NOISE{T} -> diodes, gyro, pos 

      Struct containing the noise in each sensor measurement, except for 
    some resaon I don't have magnetometer noise. This is really only used 
    for debugging anyway, so...

    I may add in plotting for this at some point...
"""
struct NOISE{N, T}
    diodes::SVector{N, T}       # Noise in each measured diode current
    gyro::SVector{3, T}         # Noise in the gyroscope measurement
    pos::SVector{3, T}          # Noise in the position estimate 
end