module Estimator



# I can just include a file for each estimator to clean up, right? 



include("estimator_config_file.jl")
# include("../system_config_file.jl")
using ..CustomStructs

using LinearAlgebra
using ForwardDiff

export estimate_vals

# Types of estimators available:
export MAG_CALIB
export DIODE_CALIB
export MEKF


include("magnetometer_calibration.jl")

# MAKE SUBMODULES...?

####################################################################
#               TRIVIAL CASE (does nothing)                        #
####################################################################
# struct TRIVIAL
#     junk
# end

function estimate_vals(sat::SATELLITE, data::TRIVIAL)
    t = data
    return sat, false 
end


# ####################################################################
# #               MAGNETOMETER CALIBRATION                           #
# ####################################################################
# # Get rid of sketchy global 

# struct MAG_CALIB 
#     mag_field_meas 
#     mag_field_pred
#     # current_meas
# end

# ######### SKETCHY 
# mag_field_hist = 0 
# curr_hist = 0 
# A = 0 
# ################


# # Currently does NOT include time-varying current-induced bias (assume we turn off everything while calibrating)
# function estimate_vals(sat::SATELLITE, data::MAG_CALIB)
#     # Sat ESTIMATE not truth 

#     # hcat(map(x -> x * I(3), [4, 5, 6, 7])...)   == kron([4,5,6,7], I(3))'
#     # Find better name for this than "currents"
#     # currents = kron(data.current_meas, I(3))' # Multiply each element by I₃

#     if isempty(size(mag_field_hist))  # Initialize 
#         # global curr_hist = data.current_meas
#         # global A = [data.mag_field_pred[1]*I(3) data.mag_field_pred[2]*I(3) data.mag_field_pred[3]*I(3) I(3) currents]

#         global mag_field_hist = data.mag_field_meas[:]    
#         global A = [(data.mag_field_pred[1]*I(3))       (data.mag_field_pred[2]*I(3))[:, 2:3]       (data.mag_field_pred[3]*I(3))[:, 3]    I(3)]
#     else
#         # global curr_hist = [curr_hist; curr_hist]
#         # new_row = [data.mag_field_pred[1]*I(3) data.mag_field_pred[2]*I(3) data.mag_field_pred[3]*I(3) I(3) currents]
        
#         global mag_field_hist = [mag_field_hist[:]; data.mag_field_meas[:]] # Store as vector
#         new_row =  [(data.mag_field_pred[1]*I(3))       (data.mag_field_pred[2]*I(3))[:, 2:3]       (data.mag_field_pred[3]*I(3))[:, 3]    I(3)]

#         global A = [A; new_row]
#     end


#     # If we have enough data to math...
#     if size(A, 1) > 30 # Needs to be overconstrained
#         params = A \ mag_field_hist 

#         # mag_calib_matrix_est = reshape(params[1:9], 3, 3)
#         # bx_est, by_est, bz_est = params[10:12]
#         # induced_scale_factors_est = reshape(params[13:end], size(sat.magnetometer.induced_scale_factors))

#         mag_calib_matrix_est, β = parameters_to_matrix_bias(params)
#         bx_est, by_est, bz_est = β[:]

#         a_est, b_est, c_est, ρ_est, λ_est, ϕ_est = extract_parameters(mag_calib_matrix_est)


#         # Check for change 
#         δscale_factors = sum(abs.(sat.magnetometer.scale_factors - [a_est, b_est, c_est]))
#         δnon_ortho = sum(abs.(sat.magnetometer.non_ortho_angles - [ρ_est, ϕ_est, λ_est]))
#         δbias = sum(abs.(sat.magnetometer.bias - [bx_est, by_est, bz_est]))
#         # δscale_coefs = sum(abs.(sat.magnetometer.induced_scale_factors - induced_scale_factors_est))

#         if (δscale_factors < 0.001) && (δnon_ortho < 0.001) && (δbias < 0.001) # && (δscale_coefs < 0.0025)
#             println("FINISHED Mag Calib!")
#             finished = true 
#         else 
#             finished = false 
#         end

#         # UPDATE SATELLITE ESTIMATES


#         updated_magnetometer_est = MAGNETOMETER([a_est, b_est, c_est], 
#                                                 [ρ_est, λ_est, ϕ_est],
#                                                 [bx_est, by_est, bz_est] )
#                                                 # induced_scale_factors_est)
#         sat.magnetometer = updated_magnetometer_est


#         # PLOT IF FINISHED??

#         return sat, finished
#     end

#     # Otherwise
#     return sat, false
# end

# function initialize(data::MAG_CALIB)
#     return data
# end

# function extract_parameters(T)
#     a = T[1,1] # Easy 

#     b = sqrt((T[2,1]^2) + (T[2,2]^2)) # (bsin)^2 + (bcos)^2 = b^2
#     ρ = atan(T[2,1], T[2,2]) # sin/cos to maintain signs

#     c = sqrt((T[3,1]^2) + (T[3,2]^2) + (T[3,3]^2))
#     ϕ = atan(T[3,2] / T[3,3])
#     λ = atan(  sign(T[3,1]) * sqrt( (T[3,1]^2) ),  
#             sign((T[3,2]^2) + (T[3,3]^2)) * sqrt( (T[3,2]^2) + (T[3,3]^2) ) ) # Not positve this portion of the sign is actually useful

#     return a, b, c, ρ, λ, ϕ
# end

# function parameters_to_matrix_bias(p)
#     # params | [9 x 1] => Lower triangular & bias vector
#     T = [p[1]   0       0;
#         p[2]   p[4]    0;
#         p[3]   p[5]    p[6]];      # Calibration matrix

#     β = p[7:9];     # Bias vector 

#     return T, β[:];
# end


####################################################################
#               DIODE CALIBRATION + MEKF                           #
####################################################################
struct DIODE_CALIB
    # albedo::ALBEDO
    # state 
    # covariance 
    # inertial_vecs 
    # ang_vel 
    # body_vecs 
    # current_meas 
    # W 
    # V

    dt 
    epc 
    num_diodes 
end

function estimate_vals(sat::SATELLITE, data::DIODE_CALIB)
    # data:
    #   state, covariance, W, V, Newtonian vectors, w, body vectors, currents, 
    #           (dt, epc, num_diodes?)

    # IF FIRST TIME -> Triad 
    # Else 
    #   Step the MEKF 
    #   Compare to SAT 
    #   Stop if low change for a bit 

    return sat, true
end

function initialize(data::DIODE_CALIB)
    # In place

    # Initialize lots of stuff 
    #   REFL data? Cell Albedos? Or pass in?
    #   ESTIMATES? Or is that done in SAT already?
    # COVARIANCE? -> Add to satellite? Or to diode data?  P, W, V
    # Initialize with TRIAD?

    return data
end




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
    return sat, true
end

function initialize(data::MEKF)
    # In place
    return data
end


end