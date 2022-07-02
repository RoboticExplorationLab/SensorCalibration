
#########
# SETUP #
#########
using Test, Infiltrator

using StaticArrays, SatelliteDynamics, EarthAlbedo 
using Distributions, LinearAlgebra, Plots, JLD2, Random
using ProgressMeter


@enum(Operation_mode, mag_cal = 1, detumble, diode_cal, mekf, chill, finished)
Base.to_index(om::Operation_mode) = Int(s)
include("../src/state_machine.jl")

function get_initial_state(; _Re = 6378136.3, detumbled = false, bias_less = false) 
    ecc = 0.0001717 + 0.00001 * randn()
    inc = 51.6426 + randn()
    Ω   = 178.1369 + randn()
    ω   = 174.7410 + randn()
    M   = 330.7918 + 100 + randn()   # +94/95 is just before sun, -40 is just before eclipse
    sma = (_Re + 421e3) / (1 + ecc)  # Apogee = semi_major * (1 + ecc)

    oe0 = [sma, ecc, inc, Ω, ω, M]   # Initial state, oscullating elements
    eci0 = sOSCtoCART(oe0, use_degrees = true) # Convert to Cartesean

    r₀ = SVector{3, Float64}(eci0[1:3])
    v₀ = SVector{3, Float64}(eci0[4:6])
    q₀ = randn(4);  q₀ = SVector{4, Float64}(q₀ / norm(q₀))
    ω₀ = (detumbled) ? SVector{3, Float64}(0.05 * randn(3)) : SVector{3, Float64}(0.5 * randn(3))
    β₀ = (bias_less) ? SVector{3, Float64}(0.0, 0.0, 0.0)  : SVector{3, Float64}(rand(Normal(0.0, deg2rad(2.0)), 3))
    
    T_orbit = orbit_period(oe0[1])
    x = STATE(r₀, v₀, q₀, ω₀, β₀)

    return x, T_orbit 
end

function get_albedo(scale = 1) 

    function load_refl(path = "data/refl.jld2", scale = 1)
        temp = load(path)
    
        refl = REFL( temp["data"][1:scale:end, 1:scale:end], temp["type"], temp["start_time"], temp["stop_time"])
    
        return refl
    end
    lat_step = 1.0 * scale
    lon_step = 1.25 * scale

    refl = load_refl("../src/data/refl.jld2", scale)  
    cell_centers_ecef = get_albedo_cell_centers(lat_step, lon_step) 
    return Simulator.ALBEDO(refl, cell_centers_ecef)
end;

function get_albedo_cell_centers(lat_step = 1, lon_step = 1.25)
    """
        Returns the cell centers for the grid covering the surface of the Earth in Cartesian ECEF, to be used in later estimations of Earth's albedo,
            by looping through each cell's LLA coordinate and converting to ECEF 

        Arguments:
        - lat_step: (Optional) The step size (in degrees) to take across the latitude domain. Defaults to 1*        | Scalar 
        - lon_step: (Optional) The step size (in degrees) to take across the longitude domain. Defaults to 1.25*    | Scalar

        Returns:
        - cells_ecef: Matrix containing [x,y,z] coordinate for each latitude, longitude point.
                        Of form [lat, lon, [x,y,z]]                                                                 | [num_lat x num_lon x 3]
    """
    alt = 0.0 # Assume all cells are on surface of earth
    num_lat = Int(round((180 - lat_step) / lat_step) + 1)
    num_lon = Int(round((360 - lon_step) / lon_step) + 1)

    lon_offset = lon_step + (360 - lon_step) / 2   # Centers at 0 (longitude: [1.25, 360] => [-179.375, 179.375])
    lat_offset = lat_step + (180 - lat_step) / 2   # Centers at 0 (latitude:  [1.00, 180] => [-89.5, 89.5])

    cells_ecef = zeros(num_lat, num_lon, 3) # Lat, Lon, [x,y,z]
    for lat = 1:num_lat 
        for lon = 1:num_lon
            geod = [(lon * lon_step - lon_offset), (lat * lat_step - lat_offset), alt]
            ecef = sGEODtoECEF(geod, use_degrees = true)

            cells_ecef[Int(lat), Int(lon), :] = ecef
        end
    end

    return cells_ecef 
end;

function rotation_noise(σ::T, dt::T) where {T}

    if σ == 0
        return (I(3))
    else 
        v   = rand(Normal(0, σ), 3)
        mag = norm(v)
        v̂   = hat(v / mag) 

        R = I(3) + (v̂) * sin(mag * dt) + (v̂^2) * (1 - cos(mag * dt))  # Equivalent to exp(v̂)

        return (R)
    end
end


##############################
##############################
##############################

@testset "Triad" begin 

    ## TRIAD - no noise
    N = 100
    es = zeros(N)
    for i = 1:N
        q = randn(4); q /= norm(q);     # Get a random rotation 
        r₁ᴵ, r₂ᴵ = randn(3), randn(3);  # Get two random vectors
        ᴮQᴵ = quat2rot(q)';
        r₁ᴮ, r₂ᴮ = ᴮQᴵ * r₁ᴵ, ᴮQᴵ * r₂ᴵ;

        q̂, ᴵQ̂ᴮ = triad(r₁ᴵ, r₂ᴵ, r₁ᴮ, r₂ᴮ)

        es[i] = norm(cayley_map(q, q̂)) < 1e-10
    end
    @test all(es .== 1)

    ## TRIAD - some noise
    N = 100
    es = zeros(N)
    for i = 1:N
        q = randn(4); q /= norm(q);     # Get a random rotation 
        r₁ᴵ, r₂ᴵ = randn(3), randn(3);  # Get two random vectors
        ᴮQᴵ = quat2rot(q)';
        η₁, η₂ = rotation_noise(deg2rad(3.0), 0.1), rotation_noise(deg2rad(3.0), 0.1)
        r₁ᴮ, r₂ᴮ = η₁ * ᴮQᴵ * r₁ᴵ, η₂ * ᴮQᴵ * r₂ᴵ;

        q̂, ᴵQ̂ᴮ = triad(r₁ᴵ, r₂ᴵ, r₁ᴮ, r₂ᴮ)

        es[i] = norm(cayley_map(q, q̂)) < 0.2
    end
    @test all(es .== 1)

    ##### Run triad #####
    N  = 100
    es = zeros(N)
    for i = 1:N
        x    = STATE()
        sat  = SATELLITE(; ideal = true)
        t    = Epoch(2021, 3, 29)
        dt   = 1.0
        alb  = get_albedo(2)

        tr, se, ecl, noise = generate_measurements(sat, alb, x, t, dt; use_albedo = false)

        q̂    = run_triad(se, sat, t, true)
        es[i] = (norm(cayley_map(x.q, q̂)) < 0.1)
    end 
    @test all(es .== 1)
end

@testset "Estimate Sun Vector" begin 

end

@testset "Correct Magnetometer" begin

end


