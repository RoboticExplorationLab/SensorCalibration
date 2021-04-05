using Pkg
using Plots
using MAT, LinearAlgebra, ForwardDiff
using JLD2
using Random, Distributions

Random.seed!(1001)

include("mekf.jl"); 
include("measurement.jl")
include("prediction.jl")
include("triad.jl")
include("rotationFunctions.jl")


@load "mekf_data.jld2"

# # load mekf_inputs
# vars = matread("TemplateCode/mekf_inputs.mat")
# rN1 = vars["rN1"][:,1]      # Pair of vectors in newtonian (inertial) frame (~noiseless?)
# rN2 = vars["rN2"][:,1]
# rB1hist = vars["rB1hist"]   # Sets of vectors in body frame (noisy)
# rB2hist = vars["rB2hist"]
# W = vars["W"]
# V = vars["V"]
# whist = vars["whist"]       # Angular Velocity (?)
# dt = vars["dt"]             # Time step

# numDiodes = 3
# Ihist = zeros(numDiodes, size(rB1hist,2))
# cVals = zeros(3,1);
# αs = zeros(3,1);
# ϵs = zeros(3,1);


i = numDiodes

μ_α = sum(αs) / i;
μ_ϵ = sum(ϵs) / i;
μ_c = sum(cVals) / i;
σ_α = deg2rad(5); # 5 degrees
σ_ϵ = deg2rad(5);
σ_c = 0.1 * μ_c;

α0 = rand(Normal(μ_α, σ_α), i)
ϵ0 = rand(Normal(μ_ϵ, σ_ϵ), i)
c0 = rand(Normal(μ_c, σ_c), i)


yhist = [rB1hist; rB2hist; Ihist]; # Measurements (bodyframe vectors)
# yhist = [rB1hist; rB2hist];

rN = hcat(rN1, rN2);

# # ADJUST THESE ##########
# W = (3.04617e-10) .* I(6)
# V = (3.04617e-4)  .* I(6) # 6 for getting rotation, i for current measurements
W = (3.04617e-10) .* I(6+3*i)
V = (3.04617e-4)  .* I(6+i) # 6 for getting rotation, i for current measurements


# Initial quaternion estimate (scalar first)
q0, R = triad(rN1,rN2,rB1hist[:,1],rB2hist[:,2]);
β0 = [0;0;0];

x0 = [q0; β0; c0; α0; ϵ0]; # Initialize with no bias, c=α=ϵ=rand, [7 x 3i]
# x0 = [q0; β0];

P0 = (10*pi/180)^2 * I((size(x0, 1) - 1)); # 10 deg. and 10 deg/sec 1-sigma uncertainty + quaternions use an extra variable, so we subtract off
#####

xhist, Phist = mekf(x0,P0,W,V,rN,whist,yhist, dt, i);


@load "mekf_truth.jld2"
# vars = matread("TemplateCode/mekf_truth.mat")
# qtrue = vars["qtrue"]
# btrue = vars["btrue"]


# Calculate error quaternions
e = zeros(size(qtrue));
for k = 1:size(qtrue,2)
    e[:,k] = qmult(qconj(qtrue[:,k]), xhist[1:4,k]); 
end






# ------ PLOTS ------ # 
### ATTITUDE 
q0Plt = plot( qtrue[4,:], linestyle = :dash, label = "q0")
q0Plt = plot!(xhist[4,:], label = false)

qiPlt = plot( qtrue[1,:], linestyle = :dash, label = "qi")
qiPlt = plot!(xhist[1,:], label = false)

qjPlt = plot( qtrue[2,:], linestyle = :dash, label = "qj")
qjPlt = plot!(xhist[2,:], label = false)

qkPlt = plot( qtrue[3,:], linestyle = :dash, label = "qk")
qkPlt = plot!(xhist[3,:], label = false)

display(plot(q0Plt, qiPlt, qjPlt, qkPlt, layout = (2,2), title = "Attitude"))

### Attitude Error 
eiPlt = plot(  (360/pi) * e[1,:], label = "Ei")
eiPlt = plot!( (360/pi) * sqrt.(Phist[1,1,:]), color = :red, label = false) 
eiPlt = plot!(-(360/pi) * sqrt.(Phist[1,1,:]), color = :red, label = false)

ejPlt = plot(  (360/pi) * e[2,:], label = "Ej")
ejPlt = plot!( (360/pi) * sqrt.(Phist[2,2,:]), color = :red, label = false) 
ejPlt = plot!(-(360/pi) * sqrt.(Phist[2,2,:]), color = :red, label = false)

ekPlt = plot(  (360/pi) * e[3,:], label = "Ek")
ekPlt = plot!( (360/pi) * sqrt.(Phist[3,3,:]), color = :red, label = false) 
ekPlt = plot!(-(360/pi) * sqrt.(Phist[3,3,:]), color = :red, label = false)

display(plot(eiPlt, ejPlt, ekPlt, layout = (3,1), title = "Attitude Error (deg)"))


### Bias Error
bxPlt = plot(  xhist[5,:] - btrue[1,:], label = "Bx")
bxPlt = plot!( 2*sqrt.(Phist[4,4,:]), color = :red, label = false)
bxPlt = plot!(-2*sqrt.(Phist[4,4,:]), color = :red, label = false)

byPlt = plot(  xhist[6,:] - btrue[2,:], label = "By")
byPlt = plot!( 2*sqrt.(Phist[5,5,:]), color = :red, label = false)
byPlt = plot!(-2*sqrt.(Phist[5,5,:]), color = :red, label = false)

bzPlt = plot(  xhist[7,:] - btrue[3,:], label = "Bz")
bzPlt = plot!( 2*sqrt.(Phist[6,6,:]), color = :red, label = false)
bzPlt = plot!(-2*sqrt.(Phist[6,6,:]), color = :red, label = false)

display(plot(bxPlt, byPlt, bzPlt, layout = (3,1), title = "Bias Error"))


xhist[(8+i):end,:] = rad2deg.(xhist[(8+i):end,:]);
αs = rad2deg.(αs)
ϵs = rad2deg.(ϵs)

### Calibration Estimates
cPlt1 = plot( xhist[8,:], color = "red", label = false)
cPlt1 = hline!([cVals[1]], color = "red", label = false, linestyle = :dash)
cPlt2 = plot(xhist[9,:], color = "blue", label = false)
cPlt2 = hline!([cVals[2]], color = "blue", label = false, linestyle = :dash)
cPlt3 = plot(xhist[10,:], color = "green", label = false)
cPlt3 = hline!([cVals[3]], color = "green", label = false, linestyle = :dash)
display(plot(cPlt1, cPlt2, cPlt3, layout = (3,1), title = "Calibration values"))

aPlt1 = plot( xhist[8+i,:], color = "red", label = false)
aPlt1 = hline!([αs[1]], color = "red", label = false, linestyle = :dash)
aPlt2 = plot(xhist[9+i,:], color = "blue", label = false)
aPlt2 = hline!([αs[2]], color = "blue", label = false, linestyle = :dash)
aPlt3 = plot(xhist[10+i,:], color = "green", label = false)
aPlt3 = hline!([αs[3]], color = "green", label = false, linestyle = :dash)
display(plot(aPlt1, aPlt2, aPlt3, layout = (3,1), title = "Azimuth Angles (α)"))

ePlt1 = plot( xhist[8+2*i,:], color = "red", label = false)
ePlt1 = hline!([ϵs[1]], color = "red", label = false, linestyle = :dash)
ePlt2 = plot(xhist[9+2*i,:], color = "blue", label = false)
ePlt2 = hline!([ϵs[2]], color = "blue", label = false, linestyle = :dash)
ePlt3 = plot(xhist[10+2*i,:], color = "green", label = false)
ePlt3 = hline!([ϵs[3]], color = "green", label = false, linestyle = :dash)
display(plot(ePlt1, ePlt2, ePlt3, layout = (3,1), title = "Elevation Angles (ϵ)"))


# display(plot(cPlt, aPlt, ePlt, layout = (3,1), title = "Calibration Estimates"))
