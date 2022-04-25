# get the virtual environment all set up
using Pkg
using Plots
# using PyPlot
cd(@__DIR__)
Pkg.activate(".")
Pkg.instantiate()


using LinearAlgebra, ForwardDiff, Attitude, StaticArrays
using SparseArrays, IterativeSolvers, Infiltrator
const FD = ForwardDiff
using SatelliteDynamics
using Printf
include("mag_field.jl")


using Random
# Random.seed!(2346)

function accel_perturbations(epc::Epoch, x::Array{<:Real} ;
             mass::Real=1.0, area_drag::Real=0.01, coef_drag::Real=2.3,
             area_srp::Real=1.0, coef_srp::Real=1.8,
             n_grav::Integer=10, m_grav::Integer=10)
    """Accelerations for spacecraft in LEO, ForwardDiff friendly"""

    # Extract position and velocity
    r = x[1:3]
    v = x[4:6]

    # Compute ECI to ECEF Transformation -> IAU2010 Theory
    PN = bias_precession_nutation(epc)
    E  = earth_rotation(epc)
    W  = polar_motion(epc)
    R  = W * E * PN

    # Compute sun and moon position
    r_sun  = sun_position(epc)
    r_moon = moon_position(epc)

    # Compute acceleration (eltype(x) makes this forward diff friendly)
    a = zeros(eltype(x), 3)

    # spherical harmonic gravity
    a += accel_gravity(x, R, n_grav, m_grav)

    # atmospheric drag
    ρ = density_harris_priester(epc,r)
    a += accel_drag([r;v],ρ,mass, area_drag, coef_drag, Array{Real, 2}(PN))

    # SRP
    nu = eclipse_conical(x, r_sun)
    a += nu*accel_srp(x, r_sun, mass, area_srp, coef_srp)

    # third body sun
    a += accel_thirdbody_sun(x, r_sun)

    # third body moon
    a += accel_thirdbody_moon(x, r_moon)

    return a
end

function measurement(x,t)
    """Measurement function for our system. Comprises magnetic field estimate along each axis"""

    r = x[1:3]
    v = x[4:6]

    t_current = epc+t*tscale

    r_sun  = sun_position(t_current)      # sun position
    η = eclipse_conical(r*dscale, r_sun)  

    currents = η * randn((numCurrents, 1))
    magFromCur = curCoef * currents;

    # Orthogonal components of the geomagnetic field
    bx, by, bz = (IGRF13(r*dscale,t_current)) 

    b_field = [bx; by; bz];
    b_field_hat = (T_matrix*b_field) + bias + magFromCur;
         
    bx_hat, by_hat, bz_hat = b_field_hat[1:3]; 

    currents /= (tscale^2) # This is just to accommodate the next line
    temp = (tscale^2)*[bx_hat; by_hat; bz_hat; 
                        bx; by; bz;
                        currents];  

    return temp[:,1];                     
end

function dynamics(x,u,t)
    """ODE for orbital motion"""
    r1 = x[1:3]*dscale
    v1 = x[4:6]*(dscale/tscale)

    # this is my debugging stuff for if one of them hits the earth
    if norm(r1)<R_EARTH
        error("Impact 1")
    end

    # the current time is (epc + t*dscale)
    t_current = epc+t*tscale

    return [v1 / (dscale/tscale);
            accel_perturbations(t_current,[r1;v1]) / (dscale/tscale^2)]
end

function rk4(f, u, x_n, h,t_n)
    """Vanilla RK4"""
    x_n = SVector{nx}(x_n)

    k1 = h*f(x_n,u,t_n)
    k2 = h*f(x_n+k1/2, u,t_n + h/2)
    k3 = h*f(x_n+k2/2, u, t_n + h/2)
    k4 = h*f(x_n+k3, u, t_n + h)


    return (x_n + (1/6)*(k1+2*k2+2*k3 + k4))

end

function generate_data(x0,T,dt,R)
    """Forward rollout of the sim, no process noise"""

    X = fill(zeros(nx),T)
    Y = fill(zeros(m),T)

    X[1] = x0

    u = SVector(0,0,0)
    for i = 1:(T-1)
        t = (i-1)*dt 
        X[i+1] = rk4(dynamics,u,X[i],dt,t)  #+ sqrt(Q)*randn(nx)
        Y[i] = measurement(X[i],t) + sqrt(R)*randn(m)
    end
    Y[T] = measurement(X[T],(T-1)*dt) +  sqrt(R)*randn(m)

    return X,Y
end

function extractParameters(T)
    a = T[1,1] # Easy 

    b = sqrt((T[2,1]^2) + (T[2,2]^2)) # (bsin)^2 + (bcos)^2 = b^2
    ρ = atan(T[2,1], T[2,2]) # sin/cos to maintain signs

    c = sqrt((T[3,1]^2) + (T[3,2]^2) + (T[3,3]^2))
    ϕ = atan(T[3,2] / T[3,3])
    λ = atan(  sign(T[3,1]) * sqrt( (T[3,1]^2) ),  
               sign((T[3,2]^2) + (T[3,3]^2)) * sqrt( (T[3,2]^2) + (T[3,3]^2) ) ) # Not positve this portion of the sign is actually useful

    return a, b, c, ρ, ϕ, λ

end


"""   SIMULATION PARAMETERS   """
T = 100 # number of knot points (~100 min= 1 orbit)

dscale, tscale = (1e6), (3600/5)
numCurrents = 25;
nx = 6                    # Size of state matrix [pos, vel]
m  = 6 + numCurrents      # Size of measurement matrix [pred, meas, current meas] 

dt = 60.0/tscale # sample time is 1 min


""" Noise """
prediction_noise =  0.0001 * (tscale^2); # Noise for for B_truth, Teslas (kg/ (amp s^2))  
measurement_noise = 0.01 * (tscale^2); # Noise for  B_meas,  Teslas (kg/ (amp s^2)) 

current_noise = 0.001;   # Noise for current magnitude vectors, Amps

R = Diagonal([prediction_noise*ones(3); measurement_noise*ones(3); current_noise*ones(numCurrents)].^2)
cholR = sqrt(R)
invcholR = inv(cholR)

""" Model Parameters """
s_a, s_b, s_c = 1 .+ 0.2*randn(3)                # Scale Factors 
bias_x, bias_y, bias_z = 10.0 * randn(3);         # Bias (μTeslas)
ρ, ϕ, λ = 15.0*randn(3);                         # Non-orthogonality angles  (in DEGREES)
ρ, ϕ, λ = deg2rad(ρ), deg2rad(ϕ), deg2rad(λ)


T_matrix = [s_a 0 0; s_b*sin(ρ) s_b*cos(ρ) 0; s_c*sin(λ) s_c*sin(ϕ)*cos(λ) s_c*cos(ϕ)*cos(λ)]
bias = [bias_x; bias_y; bias_z]

curCoef = 5*randn((3, numCurrents)) # Scales how each current affects each measurement axis 



""" Initial Conditions """
epc = Epoch(2019, 1, 1, 12, 0, 0, 0.0) # initial time for sim

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 550e3+10000*randn(), 0.01+0.01*randn(), 45.0+5*randn(), 12+2*randn(), 25+3*randn(), 0]
# oe0  = [R_EARTH + 550e3, 0.01, 45.0, 12, 25, 0]
eci0 = sOSCtoCART(oe0, use_degrees=true)

eci0[1:3] /= dscale
eci0[4:6] /= (dscale/tscale)
# println(eci0)
x0 = eci0 


# GENERATE DATA
X,Y = generate_data(x0,T,dt,R);
xVals = mat_from_vec(X)

yVals = mat_from_vec(Y)     # y_meas [3 x 1], y_pred [3 x 1], current measurements [numCurrents x 1]
y_meas = yVals[1:3, :] / (tscale^2)
y_pred = yVals[4:6, :] / (tscale^2)
currMeas = yVals[7:end, :]

# SOLVE FOR T, BIAS, and CURRENT COEFFICIENTS
d = 3
I3 = Matrix(1I, d, d)


# Clumsy method of allowing for variable number of current measurements
currRow = currMeas[1,1]*I3
for j = 2:numCurrents
    global currRow = [currRow currMeas[j,1]*I3];
end

# A = [y_pred, I, current measurements]
A = [y_pred[1,1]*I3 y_pred[2,1]*I3 y_pred[3,1]*I3 I3 currRow]


for i = 2:size(yVals)[2]
    global currRow = currMeas[1,i]*I3
    for j = 2:numCurrents
        global currRow = [currRow currMeas[j,i]*I3];
    end
    row = [y_pred[1,i]*I3 y_pred[2,i]*I3 y_pred[3,i]*I3 I3 currRow] 
    global A = [A; row]
end


y_meas_vec = reshape(y_meas, length(y_meas));

params = A \ y_meas_vec;

T_hat = reshape(params[1:9], 3, 3)
bias_hat = params[10:12]
curCoef_hat = params[13:end]
curCoef_hat = reshape(curCoef_hat, size(curCoef))

a_hat, b_hat, c_hat, ρ_hat, ϕ_hat, λ_hat = extractParameters(T_hat)
bx_hat, by_hat, bz_hat = bias_hat
ϵ_a, ϵ_b, ϵ_c, ϵ_ρ, ϵ_ϕ, ϵ_λ, ϵ_bx, ϵ_by, ϵ_bz = 
            (s_a - a_hat), (s_b - b_hat), (s_c - c_hat), (ρ -ρ_hat), (ϕ - ϕ_hat), (λ - λ_hat), (bias_x - bx_hat), (bias_y - by_hat), (bias_z - bz_hat)




println("\n\n")
println("_____________________________________________________________")
println("_____|_TRUTH____________|_PREDICTED_____|_____DIFFERENCE_____")
@printf(" a   |   %0.3f  \t| %0.3f \t|        %0.3f       \n", s_a, a_hat, ϵ_a)
@printf(" b   |   %0.3f  \t| %0.3f \t|        %0.3f       \n", s_b, b_hat, ϵ_b)
@printf(" c   |   %0.3f  \t| %0.3f \t|        %0.3f       \n", s_c, c_hat, ϵ_c)
@printf(" ρ   |   %0.3f  \t| %0.3f \t|        %0.3f       \n", ρ, ρ_hat, ϵ_ρ)
@printf(" ϕ   |   %0.3f  \t| %0.3f \t|        %0.3f       \n", ϕ, ϕ_hat, ϵ_ϕ)
@printf(" λ   |   %0.3f  \t| %0.3f \t|        %0.3f       \n", λ, λ_hat, ϵ_λ)
@printf(" bx  |   %0.3f  \t| %0.3f \t|        %0.3f       \n", bias_x, bx_hat, ϵ_bx)
@printf(" by  |   %0.3f  \t| %0.3f \t|        %0.3f       \n", bias_y, by_hat, ϵ_by)
@printf(" bz  |   %0.3f  \t| %0.3f \t|        %0.3f       \n", bias_z, bz_hat, ϵ_bz)
println("_____________________________________________________________")

println("Total T Error:      \t\t\t", sum(abs.(T_matrix - T_hat)))
println("Total Bias Error:   \t\t\t", sum(abs.(bias - bias_hat)))
println("Avg Current Coefficient Error:   \t", sum(abs.(curCoef - curCoef_hat))/length(curCoef))



magFromCur = curCoef_hat * currMeas;
y_cor = (T_hat^(-1))*((y_meas - magFromCur) .- bias)


plt = plot( y_meas[1,:], title = "Before Correction", xlabel = "Time", ylabel = "μT", color = :red, label = "Measured")
plt = plot!(y_meas[2,:], color = :blue, label = false)
plt = plot!(y_meas[3,:], color = :green, label = false)

plt = plot!(y_pred[1,:], color = :red, linestyle = :dash, label = "Predicted")
plt = plot!(y_pred[2,:], color = :blue, linestyle = :dash, label = false)
plt = plot!(y_pred[3,:], color = :green, linestyle = :dash, label = false)



pltAft = plot( y_pred[1,:], title = "After Correction", color = :red, linestyle = :dash, label = "Predicted")
pltAft = plot!(y_pred[2,:], color = :blue, linestyle = :dash, label = false)
pltAft = plot!(y_pred[3,:], color = :green, linestyle = :dash, label = false)

pltAft = plot!(y_cor[1,:], xlabel = "Time", ylabel = "μT", color = :red, label = "Adjusted")
pltAft = plot!(y_cor[2,:], color = :blue, label = false)
pltAft = plot!(y_cor[3,:], color = :green, label = false)

y_dif = abs.(y_cor - y_pred)
pltDif = plot( y_dif[1,:], color = :red, title = "Final Error")
pltDif = plot!(y_dif[2,:], color = :blue, label = false)
pltDif = plot!(y_dif[3,:], color = :green, label = false)

display(plot(pltDif))
display(plot(plt, pltAft, layout = (2,1)))