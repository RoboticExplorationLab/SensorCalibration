using Pkg
using Plots
using MAT, LinearAlgebra, ForwardDiff
using JLD2
using Random, Distributions
using BlockDiagonals           

# Random.seed!(384345)

include("mekf.jl"); 
include("sun_measurement.jl")
include("mag_measurement.jl")
include("current_measurement.jl")

include("prediction.jl")
include("triad.jl")
include("rotationFunctions.jl")

@load "mekf_data.jld2" # rB1hist rB2hist rN1 rN2 Ihist whist num_diodes dt calib_vals αs ϵs eclipse

i = num_diodes
μ_c = sum(calib_vals) / i;
σ_α = deg2rad(5); # 5 degrees
σ_ϵ = deg2rad(5);
σ_c = 0.1 * μ_c; # 10 percent

α0 = αs .+ rand(Normal(0.0, σ_α), i)
ϵ0 = ϵs .+ rand(Normal(0.0, σ_ϵ), i)
c0 = calib_vals .+ rand(Normal(0.0, σ_c), i)

yhist = [rB1hist; rB2hist; Ihist]; # Measurements (bodyframe vectors)
rN = [rN1; rN2];



#############  
estimator_params = (angle_random_walk      = 0.06,   # in deg/sqrt(hour)   
                    gyro_bias_instability  = 0.8,    # Bias instability in deg/hour
                    velocity_random_walk   = 0.014,  # in m/sec/sqrt(hour)
                    accel_bias_instability = 6)      # in microG

# Qbb = ((estimator_params[:accel_bias_instability]*(9.80665/1e6))^2)/(3600^3)
# Qgg = ((estimator_params[:gyro_bias_instability] * (pi/180)    )^2)/(3600^3)   # TODO why (seconds/hour)^3? -> Units are rad^2 /seconds^3?
# Qaa = ((estimator_params[:velocity_random_walk])^2)/3600
# Qww = ((estimator_params[:angle_random_walk]*(pi/180))^2)/3600                 # Units are rad^2 / second


Q_gyro = ((estimator_params[:gyro_bias_instability] * (pi/180)    )^2)/(3600^3)  # Units are now rad^2/seconds^3...?
σ_orient = sqrt(Q_gyro);


Q_bias = ((estimator_params[:angle_random_walk]*(pi/180))^2)/(3600)
σ_bias = sqrt(Q_bias)


Q_diode = 1e-6 # Diode Noise 

# σ_cal = 0.1 * σ_c; σ_azi = 0.1 * σ_α; σ_ele = 0.1 * σ_ϵ;
σ_cal = Q_diode; σ_azi = Q_diode; σ_ele = Q_diode;
σ_sunVec = noiseValues[:σ_η_sun]; σ_magVec = noiseValues[:σ_η_mag]; σ_curr = noiseValues[:σ_η_cur];

W = Diagonal([σ_orient * ones(3); σ_bias * ones(3); σ_cal * ones(i); σ_azi * ones(i); σ_ele * ones(i)])
V = Diagonal([σ_sunVec * ones(3); σ_magVec * ones(3); σ_curr * ones(i)])

# W = abs.(W)
# V = abs.(V)

# (Values just copied from example code)
# W = (3.04617e-10) .* I(6+3*i)
# V = (3.04617e-4)  .* I(6+i) # 6 for getting rotation, i for current measurements
#############



# Initial quaternion estimate (scalar last)
q0, R = triad(rN1[:,1],rN2[:,1],rB1hist[:,1],rB2hist[:,2]);
β0 = [0;0;0];

x0 = [q0; β0; c0; α0; ϵ0]; # Initialize with no bias, c=α=ϵ=rand, [7 x 3i]

# 10 deg, 10 deg/sec, and σ_c, σ_α, σ_ϵ 1-sigma uncertainty 
# σ_q = (10*pi/180)
# σ_β = (10*pi/180)
# p = [σ_q * ones(3); σ_β * ones(3); σ_c * ones(i); σ_α*ones(i); σ_ϵ*ones(i)].^2
# P0 = diagm(p)
P0 = (10*pi/180)^2 * I((size(x0, 1) - 1)); # 10 deg. and 10 deg/sec 1-sigma uncertainty + quaternions use an extra variable, so we subtract off

xhist, Phist = mekf(x0,P0,W,V,rN,whist,yhist, dt, i, eclipse);


@load "mekf_truth.jld2" # qtrue, btrue, pos


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
cPlt1 = plot( xhist[8,:] .- calib_vals[1], color = :blue, label = false)
cPlt1 = plot!( 2*sqrt.(Phist[7,7,:]), color = :red, label = false)
cPlt1 = plot!(-2*sqrt.(Phist[7,7,:]), color = :red, label = false)

cPlt2 = plot( xhist[9,:] .- calib_vals[2], color = :blue, label = false)
cPlt2 = plot!( 2*sqrt.(Phist[8,8,:]), color = :red, label = false)
cPlt2 = plot!(-2*sqrt.(Phist[8,8,:]), color = :red, label = false)

cPlt3 = plot( xhist[10,:] .- calib_vals[3], color = :blue, label = false)
cPlt3 = plot!( 2*sqrt.(Phist[9,9,:]), color = :red, label = false)
cPlt3 = plot!(-2*sqrt.(Phist[9,9,:]), color = :red, label = false)
display(plot(cPlt1, cPlt2, cPlt3, layout = (3,1), title = "Calibration values"))
# plot(xhist[8:(7+i),:]'); display(hline!(calib_vals, linestyle = :dash));



aPlt1 = plot( xhist[8+i,:], color = "red", label = false)
aPlt1 = hline!([αs[1]], color = "red", label = false, linestyle = :dash)
aPlt2 = plot(xhist[9+i,:], color = "blue", label = false)
aPlt2 = hline!([αs[2]], color = "blue", label = false, linestyle = :dash)
aPlt3 = plot(xhist[10+i,:], color = "green", label = false)
aPlt3 = hline!([αs[3]], color = "green", label = false, linestyle = :dash)
display(plot(aPlt1, aPlt2, aPlt3, layout = (3,1), title = "Azimuth Angles (α)"))
# plot(xhist[(8+i):(7+2*i),:]'); display(hline!(αs, linestyle = :dash));


ePlt1 = plot( xhist[8+2*i,:], color = "red", label = false)
ePlt1 = hline!([ϵs[1]], color = "red", label = false, linestyle = :dash)
ePlt2 = plot(xhist[9+2*i,:], color = "blue", label = false)
ePlt2 = hline!([ϵs[2]], color = "blue", label = false, linestyle = :dash)
ePlt3 = plot(xhist[10+2*i,:], color = "green", label = false)
ePlt3 = hline!([ϵs[3]], color = "green", label = false, linestyle = :dash)
display(plot(ePlt1, ePlt2, ePlt3, layout = (3,1), title = "Elevation Angles (ϵ)"))
# plot(xhist[(8+2*i):(7+3*i),:]'); display(hline!(ϵs, linestyle = :dash));
