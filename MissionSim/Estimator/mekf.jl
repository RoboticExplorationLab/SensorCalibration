####################################################################
#                      STANDARD MEKF                               #
####################################################################


struct MEKF 
    # state 
    # covariance 
    # inertial_vecs 
    # ang_vel 
    # body_vecs 
    # W 
    # V

    dt 
end

function estimate_vals(sat::SATELLITE, data::MEKF)
    # Initialize with end of diode_calib
    return sat, data, true
end

function initialize(data::MEKF) # Should be equal to end of diode_calib ya?
    # In place
    return data
end