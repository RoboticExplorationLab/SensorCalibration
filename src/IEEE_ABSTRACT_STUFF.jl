# Testing stuff for the IEEE Abstract 
# (right now, this means whether or not we can do mag cal and diode cal in one - AKA is it observable)

using LinearAlgebra, Distributions, ForwardDiff
include("CustomStructs.jl"); using .CustomStructs
using StaticArrays; include("quaternions.jl")


##### DYNAMICS #####

q(x) = x[1:4];   βᵧ(x) = x[5:7];
C(x) = x[8:13];  α(x) = x[14:19]; ϵ(x)  = x[20:25]
s(x) = x[26:28]; ζ(x) = x[29:31]; βₘ(x) = x[32:34]


function initial_state()
    """ x = [q, βᵧ, C, α, ϵ, s, ζ, βₘ] ∈ 𝐑⁽¹⁶⁺³ᴺ⁾ """

    q₀  = randn(4); q₀ /= norm(q₀)
    β₀  = randn(3);
    d₀  = DIODES()
    m₀  = MAGNETOMETER();
    
    return [q₀; β₀; d₀.calib_values; d₀.azi_angles; d₀.elev_angles; m₀.scale_factors; m₀.non_ortho_angles; m₀.bias]
end

# N = 500
# ω = [0.0, 0.1, -0.1]; dt = 1.0 
# xs = zeros(N, 34) 
# x = x₀
# for i = 1:N 
#     global x = dynamics(x, ω; dt = dt)
#     xs[i, :] .= x 
# end

function dynamics(x, ω⁻; dt = 1.0, σβ = 1e-3)
    ω⁺ = ω⁻ - βᵧ(x)
    nω = norm(ω⁺)

    r = ω⁺ / nω
    θ = nω * dt 
    q⁺ = (nω == 0.0) ? q(x) :
                       qmult(q(x), [cos(θ / 2); r * sin(θ / 2)])
    
    x[1:4] .= q⁺ 
    # x[5:7] += rand(Normal(0.0, σβ), 3)
    return x
end

function attitude_jacobian(q)
    G = [-q[2:4]' ; 
        q[1] * I + hat(q[2:4])]
    return G
end

function A(x, ω⁻; dt = 1.0)

    ### WITHOUT FD ###
    A₁ = zeros(33, 33)

    ω⁺ = ω⁻ - βᵧ(x)  # Corrected gyro 
    Q = exp(-hat(ω⁺) * dt)

    # Standard MEKF portion
    A₁[1:3, 1:3] .= Q
    A₁[1:3, 4:6] .= -dt * I(3)
    A₁[4:6, 4:6] .= I(3)

    # Diode portion 
    A₁[7:24, 7:24] .= I(18)

    # Magnetometer portion 
    A₁[25:33, 25:33] .= I(9)


    ### WITH FD ### 

    _dynamics(_x) = dynamics(_x, ω⁻; dt = dt, σβ = 0.0)
    A_temp = ForwardDiff.jacobian(_dynamics, x)

    A₂ = zeros(33, 33) 
    x⁺ = dynamics(x, ω⁻; dt = dt)
    q⁺ = Array(q(x⁺)) 
    q⁻ = Array(q(x ))

    A₂[1:3, 1:3] .= attitude_jacobian(q⁺)' * A_temp[1:4, 1:4] * attitude_jacobian(q⁻)   # df/dq with Attitude Jacobian 
    A₂[1:3, 4:6] .= attitude_jacobian(q⁺)' * A_temp[1:4, 5:7]
    # A₂[1:3, 7:end]  stays 0
    # A₂[4:end, 1:4]  stays 0
    A₂[4:end, 4:end] .= A_temp[5:end, 5:end]


    return A₁, A₂
end

x₀ = initial_state()
ω  = [0.0, 0.1, -0.12]
ω̃  = ω + βᵧ(x₀)
dt = 1.0 
A₁, A₂ = A(x₀, ω̃ , dt = dt)
# Differ on [1:3, 1:6], as before...


##### MEASUREMENTS #####
using SatelliteDynamics, EarthAlbedo, JLD2

include("mag_field.jl")

function mag_cal_matrix(x)
    a, b, c = s(x)
    ρ, λ, ν = ζ(x)

    T = [a         0.0              0.0;
         b*sin(ρ)  b*cos(ρ)         0.0;
         c*sin(λ)  c*sin(ν)*cos(λ)  c*cos(ν)*cos(λ)]

    return T
end

function get_albedo(scale = 1) 

    function load_refl(path = "data/refl.jld2", scale = 1)
        temp = load(path)
    
        refl = REFL( temp["data"][1:scale:end, 1:scale:end], temp["type"], temp["start_time"], temp["stop_time"])
    
        return refl
    end
    lat_step = 1.0 * scale
    lon_step = 1.25 * scale

    refl = load_refl("data/refl.jld2", scale)  
    cell_centers_ecef = get_albedo_cell_centers(lat_step, lon_step) 
    return ALBEDO(refl, cell_centers_ecef)
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

function compute_diode_albedo(data, cell_centers_ecef, surface_normal, sat_pos) 
    # @show size(data)
    Nlat, Nlon = size(data)

    diode_albedo = 0.0
    r_g = zeros(eltype(surface_normal), 3)

    for r = 1:Nlat
        for c = 1:Nlon

            if data[r,c] != 0
                r_g .= cell_centers_ecef[r, c, :] .- sat_pos
                r_g .= r_g ./ norm(r_g)  # Make unit

                cell_albedo = (data[r,c] * dot(surface_normal, r_g))

                if cell_albedo > 0.0    # Can't be negative
                    diode_albedo += cell_albedo 
                end
            end
        end
    end
    
    return diode_albedo
end

# No albedo, or noise (if doing FD)
function diode_measurement(x, alb, r, sᴵ, sᴮ; σ = 0.06, E_am₀ = 1366.9)

    Cs, αs, ϵs = C(x), α(x), ϵ(x)
    N = length(Cs)

    I  = zeros(eltype(x), N)   # Currents being generated by each photodiode
    ηI = zeros(eltype(x), N)   # Noise being applied to each photodiode
    Ĩ  = zeros(eltype(x), N)   # Noisy current measurements

    albedo_matrix = earth_albedo(SVector{3, Float64}(r), SVector{3, Float64}(sᴵ), alb.refl.data) 

    for i = 1:N 
        # Get photodiode surface normal
        surface_normal = ((cos(ϵs[i])*cos(αs[i])), (cos(ϵs[i])*sin(αs[i])), sin(ϵs[i]))   

        diode_albedo = compute_diode_albedo(albedo_matrix, alb.cell_centers_ecef, surface_normal, r) 

        # Calculate current, including noise and Earth's albedo 
        current = Cs[i] * (dot(surface_normal, sᴮ) + (diode_albedo / E_am₀) )
        η = rand(Normal(0.0, abs(σ * current)))

        # Photodiodes don't generate negative current
        I[i] = (current < 1e-8) ? 0.0 : current 
        Ĩ[i] = (current + η < 1e-8) ? 0.0 : (current + η) 
        ηI[i] = η 
    end

    return I, Ĩ, ηI
end

function get_measurements(x, t, r, alb)
    ## Generate Magnetometer measurement
    Bᴵ = IGRF13(r, t)
    Q  = quat2rot(q(x))
    Bᴮ = Q' * Bᴵ 

    Tm = mag_cal_matrix(x);
    B̃ᴮ = Tm * Bᴮ + βₘ(x) + rand(Normal(0.0, 0.15), 3);

    sᴵ  = sun_position(t);
    sᴮ  = (Q' * sᴵ)
    sᴮ /= norm(sᴮ)
    Is, Ĩ, ηI = diode_measurement(x, alb, r, sᴵ, sᴮ);

    return [Ĩ; B̃ᴮ]
end

function get_H(x, t, r, alb)
    
    ### MANUAL JACOBIAN ###
    Bᴵ = IGRF13(r, t)
    Q  = quat2rot(q(x₀))
    Bᴮ = Q' * Bᴵ 
    sᴵ  = sun_position(t);
    sᴮ  = normalize(Q' * sᴵ)
    
    H₁ = zeros(9, 33)  

    ## Magnetometer Portion
    H₁[7:9, 1:3]   .= hat(Bᴮ)      # ∂ϕ 
    H₁[7:9, 4:6]   .= zeros(3)     # ∂βᵧ
    H₁[7:9, 7:24]  .= zeros(3, 18) # ∂C,α,ϵ

    Bᴮx, Bᴮy, Bᴮz = Bᴮ
    a, b, c = s(x)
    ρ, λ, ν = ζ(x)
    H₁[7:9, 25:27] .= [Bᴮx         0                       0;    # ∂s
                       0    (Bᴮx * sin(ρ) + Bᴮy * cos(ρ))  0; 
                       0           0      (Bᴮx * sin(λ) + Bᴮy * cos(λ) * sin(ν) + Bᴮz * cos(λ) * cos(ν))  ]
    H₁[7:9, 28:30] .= [ 0   0   0;   # ∂ζ
                       (b * (Bᴮx * cos(ρ) - Bᴮy * sin(ρ)))   0   0;
                        0   (c * (Bᴮx * cos(λ) - Bᴮy * sin(λ) * sin(ν) - Bᴮz * sin(λ) * cos(ν)))   (c * (Bᴮy * cos(λ) * cos(ν) - Bᴮz * cos(λ) * sin(ν))) ]
    H₁[7:9, 31:33] .= I(3)

    ## Diode Portion
    Cs, αs, ϵs = C(x), α(x), ϵ(x)
    n = [cos.(ϵs).*cos.(αs) cos.(ϵs).*sin.(αs) sin.(ϵs)]
    ndα = [(-cos.(ϵs).*sin.(αs)) ( cos.(ϵs).*cos.(αs)) zeros(size(αs))];
    ndϵ = [(-sin.(ϵs).*cos.(αs)) (-sin.(ϵs).*sin.(αs)) cos.(ϵs)]; 
    H₁[1:6, 1:3]   .= (Cs .* n) * hat(sᴮ);           # ∂ϕ
    H₁[1:6, 4:6]   .= zeros(6, 3);                   # ∂βᵧ
    H₁[1:6, 7:12]  .= Diagonal(n * sᴮ);              # ~∂C
    H₁[1:6, 13:18] .= Diagonal(Cs .* (ndα * sᴮ));    # ∂α
    H₁[1:6, 19:24] .= Diagonal(Cs .* (ndϵ * sᴮ));    # ∂ϵ 
    H₁[1:6, 25:end] .= zeros(6, 9);                  # ∂s,ζ,βₘ


    ### FD Version ###
    _meas(_x) = get_measurements(_x, t, r, alb)
    H_temp = ForwardDiff.jacobian(_meas, x)

    H₂ = zeros(9, 33)
    q⁺ = q(x)
    H₂[:, 1:3] .=  H_temp[:, 1:4] * attitude_jacobian(q⁺)
    H₂[:, 4:end] .= H_temp[:, 5:end]

    return H₁, H₂
end


t = Epoch(2020, 1, 1)
r = [6.5e6, 1e3, -2e4]
alb = get_albedo(2);

A₁, A₂ = A(x₀, ω̃ , dt = dt);
H₁, H₂ = get_H(x₀, t, r, alb);
# A₁, H₁ = A₁[1:24, 1:24], H₁[:, 1:24];  # Just MEKF + Diode Cal 
# A₂, H₂ = A₂[1:24, 1:24], H₂[:, 1:24]; 

W0₁ = (A₁') * H₁' * H₁ * (A₁)
W0₂ = (A₂') * H₂' * H₂ * (A₂)
for i = 1:300
    A₁, A₂ = A(x₀, ω̃ , dt = dt);
    H₁, H₂ = get_H(x₀, t, r, alb);
    # A₁, H₁ = A₁[1:24, 1:24], H₁[:, 1:24];  # Just MEKF + Diode Cal 
    # A₂, H₂ = A₂[1:24, 1:24], H₂[:, 1:24];   

    global W0₁ += (A₁')^i * H₁' * H₁ * (A₁^i)
    global W0₂ += (A₂')^i * H₂' * H₂ * (A₂^i)
    
    global t += 60 * 60 * 12
    global r += [1e3, -1e3, 1e3]
end

# n  = length(x₀) - 1  # Account for q -> ϕ

# ######
A₁, H₁ = A₁[1:24, 1:24], H₁[:, 1:24];  # Just MEKF + Diode Cal 
A₂, H₂ = A₂[1:24, 1:24], H₂[:, 1:24]; 
n = 24 

# A₁, H₁ = A₁[1:6, 1:6], H₁[:, 1:6]  # Just MEKF portion
# A₂, H₂ = A₂[1:6, 1:6], H₂[:, 1:6]  
# n = 6
# ######

# l, m = size(H₁)
# O₁ = zeros(l * (n), m);
# O₂ = zeros(l * (n), m);

# for i = 1:n
#     O₁[(l * i - 8):(l * i), :] = H₁ * A₁^(i - 1);
#     O₂[(l * i - 8):(l * i), :] = H₂ * A₂^(i - 1);
# end

# @show rank(O₁)
# @show rank(O₂)


# using ControlSystems 

# B = zeros(33, 3)
# D = zeros(9, 3)
# sys = StateSpace(A₁, B, H₁, D, ControlSystems.Discrete(1.0))
# OG  = gram(sys, :o)

