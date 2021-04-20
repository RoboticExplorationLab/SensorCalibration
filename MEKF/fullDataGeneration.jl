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
include("TestingParameters/Test_MatchPaper.jl")

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

    xdot[11:13] = ((J^(-1))) * cross(-w, (J*w));

end


prob = ODEProblem(dynamics, x0, tspan, param);
sol = solve(prob, Vern7(), reltol = 1e-8, saveat = saveRate);
states = sol[:,:];

# pos = plot( sol[1,:], label = "x")
# pos = plot!(sol[2,:], label = "y")
# pos = plot!(sol[3,:], label = "z")
# display(plot(pos, title = "Position"))


# quat = plot( sol[7,:], label = "i")
# quat = plot!(sol[8,:], label = "j")
# quat = plot!(sol[9,:], label = "k")
# quat = plot!(sol[10,:], label = "Scalar")
# display(plot(quat, title = "Quaternions"))


# ang = plot( sol[11,:], label = "wx")
# ang = plot!(sol[12,:], label = "wy")
# ang = plot!(sol[13,:], label = "wz")
# display(plot(ang, title = "Angle"))



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

    ecl = eclipse_conical(pos*dscale, sN_current)

    sN_current = sN_current / norm(sN_current)
    bN_current = bN_current / norm(bN_current)


    sB = R * sN_current + η_sun_body;
    bB = R * bN_current + η_mag_body;

    Is = zeros(numDiodes)

    for i = 1:numDiodes

        n = [cos(ϵs[i])*cos(αs[i]) cos(ϵs[i])*sin(αs[i]) sin(ϵs[i])]

        Ii = cVals[i] * n * sB .+ 0; #  no albedo yet

        Is[i] = Ii[1] * ecl #   When ν = 1, current is normal, when ν = 0, no current, and in between it gets scaled (?)
    end
    if minCurrentFlag
        Is[Is .< 0] .= 0  # Photodiodes can't generate negative current
    end
    

    Y = [sB[:]; bB[:]; Is[:]];
    return Y, sN_current, bN_current, ecl
    
end


yhist = zeros((6 + numDiodes), size(states,2))


sN = zeros(3, size(states,2))
bN = zeros(3, size(states,2))
eclipse = zeros(size(states,2))
dt = saveRate

for i = 1:size(states,2)
    t = dt * (i-1);
    t_current = epc + (t * tscale)
    yhist[:,i], sN[:,i], bN[:,i], eclipse[i] = g(states[:,i], t_current)
end


whist = (states[11:13,:] + biases)/tscale  
dt = saveRate * tscale;

rB1hist = yhist[1:3,:]
rB2hist = yhist[4:6,:]
Ihist = yhist[7:end,:] + η_current

iPlt = plot( Ihist[1,:], label = false)
for i = 2:size(Ihist,1)
    global iPlt = plot!(Ihist[i,:], label = false)
end
display(plot(iPlt, title = "Measured Current"))

rPlt = plot( sN[1,:])
rPlt = plot!(sN[2,:])
rPlt = plot!(sN[3,:])
display(plot(rPlt, title = "Sun vector (N)"))

rPlt = plot( bN[1,:])
rPlt = plot!(bN[2,:])
rPlt = plot!(bN[3,:])
display(plot(rPlt, title = "Magnetic Field vector (N)"))

rN1 = sN;
rN2 = bN;

@save "mekf_data.jld2" rB1hist rB2hist rN1 rN2 Ihist whist numDiodes dt cVals αs ϵs eclipse

qtrue = states[7:10,:]
btrue = biases / tscale;

@save "mekf_truth.jld2" qtrue btrue  

println("finished")