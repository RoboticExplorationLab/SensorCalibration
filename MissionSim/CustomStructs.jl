module CustomStructs

using EarthAlbedo

export MAGNETOMETER, DIODES, SATELLITE, SENSORS, ALBEDO, GROUND_TRUTH, ESTIMATES, TRIVIAL


mutable struct MAGNETOMETER
    scale_factors::Array{<:Real, 1}          # Linear scale factors for soft iron materials (a, b, c)    |   [3,]
    non_ortho_angles::Array{<:Real, 1}       # Offset from purely orthogonal vectors (ρ, λ, ϕ)           |   [3,]
    bias::Array{<:Real, 1}                   # Constant bias on each axis (DC offset) (x₀, y₀, z₀)       |   [3,] 
    # induced_scale_factors::Array{<:Real, 2}  # Scale factors that relate current to induced mag field    |   [3 x number of diodes]
end
Base.deepcopy(s::MAGNETOMETER) = MAGNETOMETER(deepcopy(s.scale_factors), deepcopy(s.non_ortho_angles), deepcopy(s.bias)) #, deepcopy(s.induced_scale_factors))

mutable struct DIODES
    calib_values        #   [number of diodes, ]
    azi_values          #   [number of diodes, ]
    elev_values         #   [number of diodes, ]
end
Base.deepcopy(s::DIODES) = DIODES(deepcopy(s.calib_values), deepcopy(s.azi_values), deepcopy(s.elev_values))


mutable struct SATELLITE
    J #::Array{<:Real, 1}                   # Inertia Matrix of satellite  (array of 1 instead of 2 because diagonal)                              |   [3 x 3]
    magnetometer::MAGNETOMETER            # Mag calibration values                                    |   [3, 3, 3, 3 x number of diodes]
    diodes::DIODES                        # Diode calibration values 
end
Base.deepcopy(s::SATELLITE) = SATELLITE(deepcopy(s.J), deepcopy(s.magnetometer), deepcopy(s.diodes))

struct SENSORS # SHOULD THESE BE UNIT or nah?
    magnetometer  # Bᴮ
    sun           # sᴮ
    diodes
end
    
struct ALBEDO
    refl::refl_struct
    cell_centers_ecef::Array{<:Real, 3}
end

mutable struct GROUND_TRUTH # Same as sim_results?
    t_hist 
    Bᴵ_hist 
    sᴵ_hist 
    # ecl_hist
end

# mutable struct ESTIMATES
#     # B̃ᴵ_hist 
#     # s̃ᴵ_hist
#     # magnetometer::MAGNETOMETER 
#     # diodes::DIODES 
#     sat::SATELLITE
# end

struct TRIVIAL
    junk
end

end