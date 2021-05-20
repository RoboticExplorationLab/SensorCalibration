using Pkg, Plots
cd(@__DIR__)
Pkg.activate(".")  # Load in the project/manifest file stuff to get right versions
Pkg.instantiate()

using SatelliteDynamics # for mag_field.jl, Geod/ECEF conversions
using StaticArrays # for SVectors, used for speed (?)
using Random       # For generating random numbers (e.g., noise)
using JLD2         # For saving 
using LinearAlgebra 
using MATLAB       # For pulling in and working with Earth Albedo Toolbox stuff
using MAT          # For loading .mat files (could be done with matlab?)


include("mag_field.jl") # From Kevin, used for IGRF13 
include("rotationFunctions.jl") # Includes common functions for dealing with quaternions
include("satellite_configuration_file.jl") # Includes satellite parameters
# ^ Should I make this a class? Dont have to worry about clearing globals


struct refl_struct
    data
    type
    start_time
    stop_time
end

function rk4(f, params, x_n, h, t_n)
    """ 
        RK4, but with no external force, and parameters passed in instead of u.
            (Based off of code from Kevin)
    """
    x_n = SVector{num_states}(x_n)

    k1 = h*f(x_n, params, t_n)
    k2 = h*f(x_n+k1/2, params, t_n + h/2)
    k3 = h*f(x_n+k2/2, params, t_n + h/2)
    k4 = h*f(x_n+k3, params, t_n + h)

    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))
end

# NO reliance on globals (σ_β...?)
function dynamics(x, param, t) 
    """ 
        Dynamics of the system
          - x: [ [px, py, pz], [vx, vy, vz], [q⃗, q0], [wx, wy, wz], [βx, βy, βz]]     | [16,]  
          - param: Standard Gravitational Parameter, Earth radius,
                satellite orbit radius, and satellite inertial matrix                 | [4,]
          - t: Time (not currently being used, could be useful for better accel)      | Scalar

        (Note that everything is in adjusted units, not necessarily m/s)
    """
    mu, Re, r, J = param   
    xdot = zeros(size(x))

    if norm(x[1:3]) < Re # Check if satellite crashes into Earth
        error("Impact!")
    end

    # dot[px, py, pz] = [vx, vy, vz]
    xdot[1:3] = x[4:6]

    # dot[vx, vy, vz] = a⃗ = -mu r̂ / ||r||^3
    xdot[4:6] = (-mu / (norm(x[1:3])^3)) * x[1:3]

    # dot[q] = 0.5 q ⋅ [w; 0]
    w = x[11:13]
    xdot[7:10] = 0.5 * qmult(x[7:10], [w; 0])

    # dot[w] = J_inv * (-w x Jw)
    xdot[11:13] = (J^(-1)) * cross(-w, (J*w))

    # bias is just a random walk, so β = β + η
    xdot[14:16] = σ_β * randn(3)   

    return xdot
end 

# ONLY works with the default step sizes for now! SCALED Return
function get_albedo_cell_centers(lat_step = 1, lon_step = 1.25)
    """ IN DEGREES """
    alt = _Re * dscale # Assume all cells are same distance from center of earth
    num_lat = size(range(lat_step, 180, step = lat_step),1) 
    num_lon = size(range(lon_step, 360, step = lon_step),1)

    cells_ecef = zeros(num_lat, num_lon, 3) # Lat, Lon, [x,y,z]

    for lat = lat_step:lat_step:180  
        for lon = lon_step:lon_step:360 
            geod = [lon - 180.625, lat - 90.5, alt]; # Scale longitude to be [-179.375, 179.375], and latitude to be [-89.5, 89.5]
            
            lat_idx = Int(lat / lat_step)
            lon_idx = Int(lon / lon_step)
            ecef = sGEODtoECEF(geod, use_degrees = true)

            cells_ecef[lat_idx, lon_idx, :] = ecef 
        end
    end

    # return cells_ecef  # CURRENTLY UNSCALED
    return cells_ecef ./ dscale

end

# Skipping for speed
function get_diode_albedo(albedo_matrix, surface_normal, sat_pos)
    # CELL * NORMAL^T * (Dist between satellite and center of cell)
    # Create a matrix that has the [x,y,z] location of all centers of cells 
    
    cell_albedos = zeros(size(albedo_matrix))

    # SKIPPING FOR SPEED
    rows, cols = size(albedo_matrix)
    diode_albedo = 0.0
    for r = 2:3:rows
        for c = 2:3:cols
            # Is this r_g or negative r_g?
            r_g = cell_centers_ecef[r,c,:] - sat_pos # Distance from satellite to cell center
            # cell_albedos[r,c] = (albedo_matrix[r,c] * (surface_normal * r_g))[1] # Albedo contribution from cell
            diode_albedo = diode_albedo + (albedo_matrix[r,c] * (surface_normal * r_g))[1]
        end
    end

    # diode_albedo = sum(cell_albedos)
    if (diode_albedo < 0)
        diode_albedo = 0
    end
    
    return diode_albedo
end

# TODO: Add in Albedo
function measurement(x, t)  
    """
        Generates system measurements 
          - x: [ [px, py, pz], [vx, vy, vz], [q⃗, q0], [wx, wy, wz], [βx, βy, βz] ]     |  [16,]
          - t: current time since start (adjusted units)                               |  scalar
    """
    t = epc + (t * tscale)
    q_vec = x[7:9];  q_sca = x[10];

    R_B2N = I(3) + 2 * hat(q_vec) * (q_sca * I(3) + hat(q_vec)); # Equation from Kevin, quaternion -> DCM (eq 39)
    R_N2B = transpose(R_B2N)

    pos = x[1:3] * dscale
    sN = sun_position(t)         # Sun vector in Newtonian frame 
    bN = IGRF13(pos, t) # Mag vector in Newtonian frame

    ecl = eclipse_conical(pos, sN) # Determine if there is an eclipse (∈ [0, 1])

    # albedo_matrix = temp_albedo_matrix
    mat""" $albedo_matrix = albedo(-$pos, $sN, $refl); """
    # Verify that when the earth eclipses we get 0

    sN = ecl * (sN / norm(sN)) # Make unit, zero out if sun is blocked 
    bN = (bN / norm(bN))       # Make unit

    # TODO verify noise magnitude
    η_sun = get_noise_matrix(σ_η_sun, dt)
    η_mag = get_noise_matrix(σ_η_mag, dt)

    sB = η_sun * (R_N2B * sN) # (noisy) Sun vector in body frame
    bB = η_mag * (R_N2B * bN) # (noisy) Mag vector in body frame

    current_vals = zeros(num_diodes) # Currents being generated by each photodiode
    diode_albedos = zeros(num_diodes)

    posB = R_N2B * (-pos / norm(pos))

    for i = 1:num_diodes
        surface_normal = [cos(ϵs[i])*cos(αs[i]) cos(ϵs[i])*sin(αs[i]) sin(ϵs[i])]     # Photodiode surface normal 

        # (If diode isn't pointing towards EARTH then don't even bother )
        if (surface_normal * posB)[1] > 0
            diode_albedo = get_diode_albedo(albedo_matrix, surface_normal, (pos/dscale)) # earth_albedo * surface_normal' * r_g;
        else
            diode_albedo = 0
        end
        current = (calib_vals[i]  * surface_normal * sB) .+ (calib_vals[i] * diode_albedo / _E_am0) .+ σ_η_cur * randn() # Calculate current, including noise and Earth's albedo 
        current_vals[i] = current[1] * ecl  # Scale by eclipse factor 
        diode_albedos[i] = calib_vals[i] * (diode_albedo/ _E_am0)
    end

    # NOTE should I do this before or after noise?
    current_vals[current_vals .< 0] .= 0 # Photodiodes don't generate negative current

    Y =  [sB[:]; bB[:]; sN[:]; bN[:]; ecl; current_vals[:]]
    return Y, diode_albedos[:]
end

function get_noise_matrix(σ, dt)
    """
        Generates a [3 x 3] noise rotation matrix given a standard deviation 
            First generates a noise vector and then converts that into a rotation matrix
            (Note that if the standard deviation provided is 0, the identity matrix is returned)
    """
    if σ != 0.0
        η_vec = σ * randn(3) # Generate a vector 
        skew = hat(η_vec)
        norm_η = norm(η_vec)

        R = (I(3) + (skew/norm_η)*sin(norm_η*dt) + ((skew/norm_η)^2)*(1 - cos(norm_η*dt))); # Rodrigues for matrix exponential (?)
    else
        R = I(3)
    end
    return R
end

function generate_data(x0, T, dt)
    """
        Propagates the system dynamics forward for a specified number of points and updates measurment values
    """
    X = zeros(num_states, T)
    Y = zeros(num_meas, T)
    X[:,1] = x0

    diode_albedos = zeros(num_diodes, T)

    for i = 1:(T-1)
        t = (i - 1) * dt  
        X[:, i + 1] = rk4(dynamics, dynamics_params, X[:,i], dt, t)
        Y[:, i], diode_albedos[:,i] = measurement(X[:,i], t)
        if i % 50 == 0
            @show i
        end
    end
    Y[:,T], diode_albedos[:,T] = measurement(X[:,T], (T-1)*dt)

    X[11:13,:] = X[11:13,:] .+ X[14:16,:] # Add in bias to the omega vectors (Note that this is done here so as to not interfere with the dynamics)
    return X, Y, diode_albedos
end

# Load in reflectivity data as a struct 
cell_centers_ecef = get_albedo_cell_centers()
mat""" addpath TEMP """
refl_dict = matread("refl.mat")
# temp_albedo_matrix = matread("albedo_matrix.mat")
# temp_albedo_matrix = temp_albedo_matrix["albedo_matrix"]
refl = refl_struct(refl_dict["data"], refl_dict["type"], refl_dict["start_time"], refl_dict["stop_time"])

###### PROPAGATE DYNAMICS #####
states, meas, albedos = generate_data(x0, T, dt)
println("Finished propagating dynamics")


##### GENERATE MEASUREMENTS #####
# History of respective values 
sB_hist = zeros(3, T)
bB_hist = zeros(3, T)
I_hist  = zeros(num_diodes, T)
sN_hist = zeros(3, T) 
bN_hist = zeros(3, T) 
ecl_hist= zeros(T)    

for i = 1:T 
    sB_hist[:,i] = meas[1:3,i] 
    bB_hist[:,i] = meas[4:6,i] 
    sN_hist[:,i] = meas[7:9,i] 
    bN_hist[:,i] = meas[10:12,i]  
    ecl_hist[i]  = meas[13,i]
    I_hist[:,i] = meas[14:(13+num_diodes),i]
end
println("Finished splitting measurements")


##### SAVE DATA #####
whist = states[11:13,:] / tscale 
biases = states[14:16,:] / tscale
dt = dt * tscale
noiseValues = (σ_η_sun = σ_η_sun, σ_η_mag = σ_η_mag, σ_η_cur = σ_η_cur)

rB1hist = sB_hist; rB2hist = bB_hist; rN1 = sN_hist; rN2 = bN_hist; eclipse = ecl_hist; Ihist = I_hist

@save "mekf_data.jld2" rB1hist rB2hist rN1 rN2 Ihist whist num_diodes dt calib_vals αs ϵs eclipse noiseValues

qtrue = states[7:10,:]
btrue = biases
pos = states[1:3,:] * dscale

@save "mekf_truth.jld2" qtrue btrue pos
println("Finished Saving")

##### PLOTS #####
display(plot(states[1:3,:]', title = "Positions"))
display(plot(sN_hist', title = "Sun Vector (N)"))
display(plot(bN_hist', title = "Mag Vector (N)"))
display(plot(I_hist', title = "Diode Currents (A)"))

display(plot(states[14:16,:]', title = "Biases"))
display(plot(states[11:13,:]', title = "W"))

display(plot(albedos', title = "Albedos"))


# using MAT
# println("UNNORMALIZED AND UNSAVED")
# matwrite("unscaled_sat_sun_vectors.mat", Dict("sun" => sN_hist, "sat" => pos); compress = true)


