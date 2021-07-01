module CustomStructs

using EarthAlbedo, SatelliteDynamics

export MAGNETOMETER, DIODES, SATELLITE, SENSORS, ALBEDO, GROUND_TRUTH, ESTIMATES, TRIVIAL

# Redo all in alphabetical order...? Or can you name it when setting up like DIODES(calib = 4, azi = 3)

mutable struct MAGNETOMETER
    scale_factors::Array{<:Real, 1}          # Linear scale factors for soft iron materials (a, b, c)    |   [3,]
    non_ortho_angles::Array{<:Real, 1}       # Offset from purely orthogonal vectors (ρ, λ, ϕ)           |   [3,]
    bias::Array{<:Real, 1}                   # Constant bias on each axis (DC offset) (x₀, y₀, z₀)       |   [3,] 
end
Base.deepcopy(s::MAGNETOMETER) = MAGNETOMETER(deepcopy(s.scale_factors), deepcopy(s.non_ortho_angles), deepcopy(s.bias)) #, deepcopy(s.induced_scale_factors))

mutable struct DIODES
    calib_values::Array{Float32, 1}      #   [number of diodes, ]
    azi_angles::Array{Float32, 1}         #   [number of diodes, ]
    elev_angles::Array{Float32, 1}       #   [number of diodes, ]
end
Base.deepcopy(s::DIODES) = DIODES(deepcopy(s.calib_values), deepcopy(s.azi_angles), deepcopy(s.elev_angles))


mutable struct SATELLITE
    J #::Array{<:Real, 1}                   # Inertia Matrix of satellite  (array of 1 instead of 2 because diagonal)                              |   [3 x 3]
    magnetometer::MAGNETOMETER            # Mag calibration values                                    |   [3, 3, 3, 3 x number of diodes]
    diodes::DIODES                        # Diode calibration values 
end
Base.deepcopy(s::SATELLITE) = SATELLITE(deepcopy(s.J), deepcopy(s.magnetometer), deepcopy(s.diodes))

mutable struct SENSORS # SHOULD THESE BE UNIT or nah?
    magnetometer::Array{Float32, 1}  # Bᴮ
    sun::Array{Float32, 1}           # sᴮ
    diodes::Array{Float32, 1}
    gyro::Array{Float32, 1}
    gps::Array{Float32, 1}           # Position needed for albedo. Perfect for now
end
    
struct ALBEDO
    refl::refl_struct
    cell_centers_ecef::Array{<:Real, 3}
end

mutable struct GROUND_TRUTH # Same as sim_results?
    t_hist::Epoch
    Bᴵ_hist::Array{Float32, 1}
    sᴵ_hist::Array{Float32, 1}
end


struct TRIVIAL
    junk
end

end