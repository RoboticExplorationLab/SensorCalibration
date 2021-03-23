using Pkg
using Plots
using MAT, LinearAlgebra, ForwardDiff

include("mekf.jl"); 
include("measurement.jl")
include("prediction.jl")
include("triad.jl")
include("rotationFunctions.jl")

i = 2
α0 = 5 * randn(i, 1) * pi / 180 
ϵ0 = 5 * randn(i, 1) * pi / 180 
c0 = 5 * randn(i, 1) 
wc, wα, wϵ = 1, 1, 1

# load mekf_inputs
vars = matread("mekf_inputs.mat")
rN1 = vars["rN1"][:,1]      # Pair of vectors in newtonian (inertial) frame (~noiseless?)
rN2 = vars["rN2"][:,1]
rB1hist = vars["rB1hist"]   # Sets of vectors in body frame (noisy)
rB2hist = vars["rB2hist"]
W = vars["W"]
WVals =  diag(W)
WVals = [WVals; ones(i,1)*wc; ones(i,1)*wα; ones(i,1)*wϵ]
W = Diagonal(WVals[:,1])


V = vars["V"]
whist = vars["whist"]
dt = vars["dt"]             # Time step

rN = hcat(rN1, rN2);
yhist = [rB1hist; rB2hist]; # Measurements (bodyframe vectors)

# Initial quaternion estimate (scalar first)
q0, R = triad(rN1,rN2,rB1hist[:,1],rB2hist[:,2]);
β0 = [0;0;0];
ang0 = [α0; ϵ0];

x0 = [q0; β0; c0; α0; ϵ0]; # Initialize with no bias, α=ϵ=0


# Correct dimensions for P?
P0 = (10*pi/180)^2*I((size(x0, 1) - 1)); # 10 deg. and 10 deg/sec 1-sigma uncertainty + quaternions use an extra variable, so we subtract off

xhist, Phist = mekf(x0,P0,W,V,rN,whist,yhist, dt, i);

# load mekf_truth
vars = matread("mekf_truth.mat")
btrue = vars["btrue"]
qtrue = vars["qtrue"]

# qtrue is scalar last, so we need to rearrange it  ##########
# temp = qtrue 
# qtrue[1,:] = temp[4,:]
# qtrue[2:4,:] = temp[1:3,:]

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
eiPlt = plot!( (360/pi) * sqrt.(Phist[1,1,:]), color = :red, label = "Variance bound?") # Should I squeeze Phist?
eiPlt = plot!(-(360/pi) * sqrt.(Phist[1,1,:]), color = :red, label = false)

ejPlt = plot(  (360/pi) * e[2,:], label = "Ej")
ejPlt = plot!( (360/pi) * sqrt.(Phist[2,2,:]), color = :red, label = "Variance bound?") # Should I squeeze Phist?
ejPlt = plot!(-(360/pi) * sqrt.(Phist[2,2,:]), color = :red, label = false)

ekPlt = plot(  (360/pi) * e[3,:], label = "Ek")
ekPlt = plot!( (360/pi) * sqrt.(Phist[3,3,:]), color = :red, label = "Variance bound?") # Should I squeeze Phist?
ekPlt = plot!(-(360/pi) * sqrt.(Phist[3,3,:]), color = :red, label = false)

display(plot(eiPlt, ejPlt, ekPlt, layout = (3,1), title = "Attitude Error (deg)"))


### Bias Error
bxPlt = plot(  xhist[5,:] - btrue[1,:], label = "Bx")
bxPlt = plot!( 2*sqrt.(Phist[4,4,:]), color = :red, label = "Variance?")
bxPlt = plot!(-2*sqrt.(Phist[4,4,:]), color = :red, label = false)

byPlt = plot(  xhist[6,:] - btrue[2,:], label = "By")
byPlt = plot!( 2*sqrt.(Phist[5,5,:]), color = :red, label = "Variance?")
byPlt = plot!(-2*sqrt.(Phist[5,5,:]), color = :red, label = false)

bzPlt = plot(  xhist[7,:] - btrue[3,:], label = "By")
bzPlt = plot!( 2*sqrt.(Phist[6,6,:]), color = :red, label = "Variance?")
bzPlt = plot!(-2*sqrt.(Phist[6,6,:]), color = :red, label = false)

display(plot(bxPlt, byPlt, bzPlt, layout = (3,1), title = "Bias Error"))