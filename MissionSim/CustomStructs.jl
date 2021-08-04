module CustomStructs

using EarthAlbedo, SatelliteDynamics

export MAGNETOMETER, DIODES, SATELLITE, SENSORS, ALBEDO, GROUND_TRUTH, ESTIMATES, TRIVIAL, FLAGS, NOISE

# Redo all in alphabetical order...? Or can you name it when setting up like DIODES(calib = 4, azi = 3)

mutable struct MAGNETOMETER
    scale_factors::Array{<:Real, 1}          # Linear scale factors for soft iron materials (a, b, c)    |   [3,]
    non_ortho_angles::Array{<:Real, 1}       # Offset from purely orthogonal vectors (ρ, λ, ϕ)           |   [3,]
    bias::Array{<:Real, 1}                   # Constant bias on each axis (DC offset) (x₀, y₀, z₀)       |   [3,] 
end
Base.deepcopy(s::MAGNETOMETER) = MAGNETOMETER(deepcopy(s.scale_factors), deepcopy(s.non_ortho_angles), deepcopy(s.bias)) #, deepcopy(s.induced_scale_factors))

mutable struct DIODES
    calib_values::Array{Float64, 1}      #   [number of diodes, ]
    azi_angles::Array{Float64, 1}         #   [number of diodes, ]
    elev_angles::Array{Float64, 1}       #   [number of diodes, ]
end
Base.deepcopy(s::DIODES) = DIODES(deepcopy(s.calib_values), deepcopy(s.azi_angles), deepcopy(s.elev_angles))


mutable struct SATELLITE
    J #::Array{Float64, 2}                   # Inertia Matrix of satellite  (array of 1 instead of 2 because diagonal)    |   [3 x 3]
    magnetometer::MAGNETOMETER            # Mag calibration values      |   [3, 3, 3, 3 x number of diodes]
    diodes::DIODES                        # Diode calibration values 
    state #::Array{<:Real, 1}               # Satellite State             | [q⃗ q₀ β⃗]
end
Base.deepcopy(s::SATELLITE) = SATELLITE(deepcopy(s.J), deepcopy(s.magnetometer), deepcopy(s.diodes), deepcopy(s.state))

mutable struct SENSORS # SHOULD THESE BE UNIT or nah?
    magnetometer::Array{Float64, 1}  # Bᴮ
    diodes::Array{Float64, 1}
    gyro::Array{Float64, 1}
    gps::Array{Float64, 1}           # Position needed for albedo. Perfect for now
end

struct NOISE 
    magnetometer::Array{Float64, 1} 
    diodes::Array{Float64, 1}
    gyro::Array{Float64, 1}
    gps::Array{Float64, 1}
    sᴮ_rot::Array{Float64, 2}
    Bᴮ_rot::Array{Float64, 2}
end
    
struct ALBEDO
    refl::refl_struct
    cell_centers_ecef::Array{<:Real, 3}
end

mutable struct GROUND_TRUTH # Same as sim_results?
    t_hist::Epoch
    Bᴵ_hist::Array{Float64, 1}
    sᴵ_hist::Array{Float64, 1}
    ŝᴮ_hist::Array{Float64, 1}
    Bᴮ_hist::Array{Float64, 1}
end


struct TRIVIAL
    junk
end


mutable struct FLAGS
    in_sun::Bool 
    magnetometer_calibrated::Bool
    diodes_calibrated::Bool
    detumbling::Bool
    calibrating::Bool
end


end