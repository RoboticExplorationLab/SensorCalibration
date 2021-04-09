using Pkg
using Plots
cd(@__DIR__)
Pkg.activate(".")
Pkg.instantiate()

using DifferentialEquations
using LinearAlgebra
using Random, Distributions
using Attitude, SatelliteDynamics
using JLD2
include("mag_field.jl")

# Random.seed!(2534)

dscale = 1#1e6
tscale = 1#3600/5

# Constants
_Re = 6378136.3 / dscale            # Radius of the earth 
_mu = (tscale^2)*(3.9860044188)*(10^(14)) / (dscale ^ 3);  # (m^3/s^2) Standard Gravitational Parameter
_r0 = (550000 / dscale) + _Re       # Distance from two center of masses
_a = _r0                            # Semi-major axis of (elliptical) orbit (circular orbit if a = r0)
_v0 = sqrt(_mu*( (2/_r0) - (1/_a))) 

# Initial Quaternion
r = [1; 1; 1/2]; r = r / norm(r);
θ = pi/4;
q0 =  [r * sin(θ/2); cos(θ/2)];

# Initial position, velocity, and angular velocity 
p0 = [_r0, 0, 0];
v0 = [0, _v0, 0];
# w0 = [0.0001, 0.0002, 0.0003] * tscale;
w0 = [0.0002, 0.0004, 0.0006] * tscale;

x0 = [p0; v0; q0; w0];
param = [_mu, _Re, _r0]; # Parameters


# Simulator Parameters
orbitTime = ((2*pi*_r0)/(_v0));
day = 24.0 * 60 * 60 / tscale;
tspan = (0, 3 * orbitTime) 
# tspan = (0, 60*150 / tscale)
saveRate = 1 / tscale;    


# State Propagation
function qmult(q, p)
    # Performs quaternion multiplication (Hamilton product) 
    #   (Quaternions are scalar last)

    q = q[:]
    p = p[:]
    
    q0, p0 = q[4], p[4]
    q⃗ = q[1:3]
    p⃗ = p[1:3]

    r0 = q0 * p0 - dot(q⃗, p⃗)
    r⃗ = (q0 * p⃗) + (p0 * q⃗) + cross(q⃗, p⃗)

    r = [r⃗; r0]
    return r
end

function dynamics(xdot, x, p, t)
    # x = [px py pz vx vy vz quat(scalar last) wx wy wz] + bias terms

    mu, Re, r = p
    xdot[1:3] = x[4:6]; 

    normr = norm(x[1:3])
    a = (-mu / (normr^3)) * x[1:3]
    xdot[4:6] = a 

    q = x[7:10]

    w = x[11:13]

    xdot[7:10] = 0.5 * qmult(q, [w; 0])

    J = [1 0 0; 0 2 0; 0 0 3];   # # # # # # # PLACEHOLDER
    xdot[11:13] = ((J^(-1))) * cross(-w, (J*w));

end


prob = ODEProblem(dynamics, x0, tspan, param);
sol = solve(prob, Vern7(), reltol = 1e-8, saveat = saveRate);

epc = Epoch(2019, 1, 1, 12, 0, 0, 0.0) # initial time for sim

pos = plot( sol[1,:], label = "x")
pos = plot!(sol[2,:], label = "y")
pos = plot!(sol[3,:], label = "z")
display(plot(pos, title = "Position"))


quat = plot( sol[7,:], label = "i")
quat = plot!(sol[8,:], label = "j")
quat = plot!(sol[9,:], label = "k")
quat = plot!(sol[10,:], label = "Scalar")
display(plot(quat, title = "Quaternions"))


ang = plot( sol[11,:], label = "wx")
ang = plot!(sol[12,:], label = "wy")
ang = plot!(sol[13,:], label = "wz")
display(plot(ang, title = "Angle"))


states = sol[:,:];


# Measurement
function hat(x)
    x = x[:];
    h = [0     -x[3]  x[2];
         x[3]    0   -x[1];
         -x[2]  x[1]    0];
end

function quat2rot(quat)
    # Scalar-last quaternions
    
    q = quat[:];
    
    s = norm(q)^(-2);
    qr = q[4];
    qi = q[1];
    qj = q[2];
    qk = q[3];
    
    R = [(1 - 2*s*(qj^2 + qk^2)) (2*s*(qi*qj - qk*qr)) (2*s*(qi*qk + qj*qr));
        (2*s*(qi*qj + qk*qr))   (1-2*s*(qi^2+qk^2))   (2*s*(qj*qk - qi*qr));
        (2*s*(qi*qk - qj*qr))   (2*s*(qj*qk+qi*qr))   (1-2*s*(qi^2 + qj^2))];
    
end

function g(x, t) # Epc has already been added to t
    pos = x[1:3];
    q = x[7:10];

    qVec = q[1:3];
    qSca = q[4];
    R = I(3) + 2 * hat(qVec) * (qSca * I(3) + hat(qVec)); # R N->B   
    R = transpose(R); # R B->N

    sN_current = sun_position(t)   
    bN_current = IGRF13(pos*dscale, t)
    sN_current = sN_current / norm(sN_current)
    bN_current = bN_current / norm(bN_current)


    sB = R * sN_current + 0.005 * randn(3);
    bB = R * bN_current + 0.01 * randn(3);

    Is = zeros(numDiodes)

    for i = 1:numDiodes

        n = [cos(ϵs[i])*cos(αs[i]) cos(ϵs[i])*sin(αs[i]) sin(ϵs[i])]

        Ii = cVals[i] * n * sB .+ 0 .+ 0.01*randn(1); # ση = 0.1, no albedo yet

        Is[i] = Ii[1]
    end

    Y = [sB[:]; bB[:]; Is[:]];
    return Y, sN_current, bN_current
    
end


numDiodes = 5;
cVals = 1 .+ 0.2 * randn(numDiodes);

ϵs = (pi/4) .+ 0.2 * randn(numDiodes);
#   Elevation ∈ [0, π/2]
idxs = (ϵs .> (pi/2));
ϵs[idxs] .= (pi/2);
idxs = (ϵs .< 0);
ϵs[idxs] .= 0


# Azimuth ∈ [0, 2π)
αs = (pi/2) .+ 0.3 * randn(numDiodes);
idxs = (αs .>= (2*pi));
αs[idxs] .= (2*pi - 0.1);
idxs = (αs .< 0);
αs[idxs] .= 0


yhist = zeros((6 + numDiodes), size(states,2))


# sN_c = [  0.6692628932588212;  # Manchesters
#         -0.681848889145605;
#         0.2952444276827865];

# bN_c = [  0.26270092455049937;
#         -0.7082385474829961;
#         -0.6552758076562027];

# sN = [0.19126880387025275; # Averaged
#      -0.9005430099234119; 
#      -0.3904331780119909];

# bN = [0.28322329459216955;
#      -0.04528518125275634;  
#       0.9574935528139525];


sN = zeros(3, size(states,2))
bN = zeros(3, size(states,2))
dt = saveRate

for i = 1:size(states,2)
    t = dt * (i-1);
    t_current = epc + (t * tscale)
    yhist[:,i], sN[:,i], bN[:,i] = g(states[:,i], t_current)
end


β0 = w0 * 0.1;
b1 = rand(Normal(β0[1], 0.1 * β0[1]), size(yhist,2))'
b2 = rand(Normal(β0[2], 0.1 * β0[2]), size(yhist,2))'
b3 = rand(Normal(β0[3], 0.1 * β0[3]), size(yhist,2))'

biases = [b1; b2; b3];
whist = (states[11:13,:]+ biases)/tscale   # Not sure if i added bias in correctly. Should I just do it here...?
biases = biases / tscale
dt = saveRate * tscale;

rB1hist = yhist[1:3,:]
rB2hist = yhist[4:6,:]
Ihist = yhist[7:end,:]

iPlt = plot( yhist[7,:])
iPlt = plot!(yhist[8,:])
iPlt = plot!(yhist[9,:])
display(plot(iPlt, title = "Measured Current"))

rPlt = plot( sN[1,:])
rPlt = plot!(sN[2,:])
rPlt = plot!(sN[3,:])
display(plot(rPlt, title = "Sun vector (N)"))

rPlt = plot( bN[1,:])
rPlt = plot!(bN[2,:])
rPlt = plot!(bN[3,:])
display(plot(rPlt, title = "Other vector (N)"))

rN1 = sN;
rN2 = bN;

@save "mekf_data.jld2" rB1hist rB2hist rN1 rN2 Ihist whist numDiodes dt cVals αs ϵs 

qtrue = states[7:10,:]
btrue = biases;

@save "mekf_truth.jld2" qtrue btrue

println("finished")