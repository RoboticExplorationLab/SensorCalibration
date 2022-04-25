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


using Random, Distributions
# Random.seed!(8236)

function hat(v)
    # Forms a skew symmetric matrix from a vector 

    v = v[:];
    vhat = [  0   -v[3]  v[2]; 
             v[3]   0   -v[1];
            -v[2]  v[1]   0];

    return vhat;
end

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


    # Orthogonal components of the geomagnetic field
    bx, by, bz = (IGRF13(r*dscale,t_current)) 

    b_field = [bx; by; bz];

    k = randn(3); k = k / norm(k) # make a unit vector 
    K = hat(k)

    θ = rand(Uniform(0, 2*pi))[1]
    R = I(3) + sin(θ) * K + (1 - cos(θ)) * K^2

    # R = I(3)

    b_fieldᴮ = R * b_field
    b_field_hat = (T_matrix*b_fieldᴮ) + bias;
         
    bx_hat, by_hat, bz_hat = b_field_hat[1:3]; 

    temp = (tscale^2)*[bx_hat; by_hat; bz_hat; 
                        bx; by; bz;
                        b_fieldᴮ[:]];  

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

function params_to_matrix(p) # Make a function that splits params into Bias and T ??
    # params | [6 x 1] => Lower triangular 
    T = [p[1]   0       0;
         p[2]   p[4]    0;
         p[3]   p[5]    p[6]];

    return T
end

function f(bm, p)
    # Reshape p -> T, b, s 
    T_hat =  params_to_matrix(p[1:6])
    bias_hat = p[7:9]

    # B_meas = TB + b + Sum(s*I_meas)
    # -> B = T^-1(B_meas - b - sum)
    B = (T_hat^(-1))*(bm - bias_hat)
    B_squared = (B[1]^2) + (B[2]^2) + (B[3]^2)
    # return (B.^2)
    return B_squared
end

function f2(be, cm, p)
    # Reshape p -> T, b, s 
    T_hat = reshape(p[1:9], 3, 3)
    bias_hat = p[10:12]
    curCoef_hat = p[13:end]
    curCoef_hat = reshape(curCoef_hat, size(curCoef))

    # B_meas = TB + b + Sum(s*I_meas)
    # -> B = T^-1(B_meas - b - sum)
    magFromCur = curCoef_hat * cm;
    B = T_hat*be + bias_hat + magFromCur
    return (B.^2)
end

function residual(x)
    """ residual vector for Gauss-Newton. rᵀr = MAP cost function
            Note that x = parameters 
            Meas = [(B_meas, B_pred) x T]
            Loss Function: 
               J = 0.5*(B^2 - f(B,I,x))^T(B^2 - f(B,I,x))
    """
    r = zeros(eltype(x), (T))
    
    for i = 1:T 
        B_meas = y_meas[:,i]
        # B_exp = y_pred[:, i] 
        B_exp_squared = (y_pred[1,i]^2) + (y_pred[2,i]^2) + (y_pred[3,i]^2) # Should be [M x 1] -> Should M > 1?


        # J = ((B_exp_squared) - f(B_meas, x)) 
        # J = 0.5 * (B_exp_squared - f(B_meas, x))' * (B_exp_squared - f(B_meas, x))

        J = f(B_meas, x) - B_exp_squared
        r[i] = invcholQ * J
    end

    return reshape(r, length(r))
end

function sparse_jacobian!(J,x)
    """Modify sparse jacobian in place"""
    #NOTE: this is just a fancy version of FD.jacobian(residual,x)
    u = 0.0
    t = 0.0
    _C_fx(x) = measurement(x,t)
    for i = 1:T


        if i < T
            k = idx_x[i]
            xt = @view x[k]
            t = (i-1)*dt
            A_fx(x) = rk4(dynamics, u, x, dt,t)
            J[k,k] = -invcholQ*FD.jacobian(A_fx,xt)
            J[k,k .+ nx] = invcholQ
            J[idx_y[i],k] = -invcholR*FD.jacobian(_C_fx,xt)
        else
            k = idx_x[i]
            xt = @view x[k]
            t = (i-1)*dt
            J[idx_y[i],k] = -invcholR*FD.jacobian(_C_fx,xt)
        end
    end
    return nothing
end

function gauss_newton(x0)
    """Gauss-Newton for batch estimation"""

    # copy initial guess
    x = copy(x0)

    # create sparse jacobian
    # J = spzeros(nx*(T-1) + m*(T),nx*T)

    Ds = 0.0
    v = zeros(length(x))

    # run Gauss-Newton for 100 iterations max
    for i = 1:50

        # ∂r/∂x
        # sparse_jacobian!(J,x)
        # this is the same as:
        # _C_fx(x) = residual(x, yVals)
        J = FD.jacobian(residual,x)

        # calculate residual at x
        r = residual(x)

        # solve for Gauss-Newton step (direct, or indirect)
        v = -J\r
        # lsqr!(v,-J,r)

        # calculate current cost
        S_k = dot(r,r)
        # @show S_k

        # step size (learning rate)
        α = 1.0

        # run a simple line search
        for ii = 1:25
            x_new = x + α*v
            S_new= norm(residual(x_new))^2

            # this could be updated for strong frank-wolfe conditions
            if S_new < S_k
                x = copy(x_new)
                Ds = S_k - S_new
                # @show ii
                break
            else
                α /= 2
            end
            if ii == 25
                @warn "line search failed"
                Ds = 0

            end
        end

        # depending on problems caling, termination criteria should be updated
        if Ds < 1e-5
            break
        end

        # ----------------------------output stuff-----------------------------
        # if rem((i-1),4)==0
        #     println("iter      α           S          dS")
        # end
        # S_display = round(S_k,sigdigits = 3)
        # dS_display = round(Ds,sigdigits = 3)
        # alpha_display = round(α,sigdigits = 3)
        # println("$i         $alpha_display   $S_display    $dS_display")

    end
    return x
end




scale_flag = false
error_flag = false
"""   SIMULATION PARAMETERS   """
T = 100 # number of knot points (~100 = 1 orbit)

dscale, tscale = 1, 1 #(1e6), (3600/5)
nx = 6                    # Size of state matrix [pos, vel]
m  = 9                    # Size of measurement matrix [pred, meas] 

dt = 60.0/tscale  # sample time is 1 min


""" Noise """
prediction_noise =  0.0001 * (tscale^2); # Noise for for B_truth, Teslas (kg/ (amp s^2))  
measurement_noise = 0.01 * (tscale^2); # Noise for  B_meas,  Teslas (kg/ (amp s^2)) 

R = Diagonal([prediction_noise*ones(3); measurement_noise*ones(3); measurement_noise*ones(3)].^2)
cholR = sqrt(R)
invcholR = inv(cholR)

Q = (1)*Diagonal(@SVector ones(3))   # Arbitrarily chosen (and not used...?)
cholQ = sqrt(Q)
invcholQ = 1 # inv(cholQ) -> Doesn't seem to make a difference

""" Model Parameters """
s_a, s_b, s_c = 1 .+ 0.2*randn(3)                # Scale Factors 
bias_x, bias_y, bias_z = 10.0 * randn(3);         # Bias (Teslas)
ρ, ϕ, λ = 15.0*randn(3);                         # Non-orthogonality angles  (in DEGREES)
ρ, ϕ, λ = deg2rad(ρ), deg2rad(ϕ), deg2rad(λ)


T_matrix = [s_a 0 0; s_b*sin(ρ) s_b*cos(ρ) 0; s_c*sin(λ) s_c*sin(ϕ)*cos(λ) s_c*cos(ϕ)*cos(λ)]
bias = [bias_x; bias_y; bias_z]

true_params = [s_a; 0; 0; 
            s_b*sin(ρ); s_b*cos(ρ); 0; 
            s_c*sin(λ); s_c*sin(ϕ)*cos(λ); s_c*cos(ϕ)*cos(λ); bias[:]]


""" Initial Conditions """
epc = Epoch(2019, 1, 1, 12, 0, 0, 0.0) # initial time for sim

# Declare initial state in terms of osculating orbital elements
oe0  = [R_EARTH + 550e3+10000*randn(), 0.01+0.01*randn(), 45.0+5*randn(), 12+2*randn(), 25+3*randn(), 0]
# oe0  = [R_EARTH + 550e3, 0.01, 45.0, 12, 25, 0]
eci0 = sOSCtoCART(oe0, use_degrees=true)

eci0[1:3] /= dscale
eci0[4:6] /= (dscale/tscale)
x0 = eci0 


# GENERATE DATA
X,Y = generate_data(x0,T,dt,R);
xVals = mat_from_vec(X)

yVals = mat_from_vec(Y)     # y_meas [3 x 1], y_pred [3 x 1]
y_meas = yVals[1:3, :] / (tscale^2)
y_pred = yVals[4:6, :] / (tscale^2)
yᴮ_true = yVals[7:9, :] / (tscale^2)

# SOLVE FOR T, BIAS
#   x 0 0 0 0 0 1 0 0
#   0 x 0 y 0 0 0 1 0           => Basically just trimming off the an extra column in y and two in z (from the y * I₃
#   0 0 x 0 y z 0 0 1

A = [(y_pred[1,1]*I(3))   (y_pred[2,1]*I(3))[:, 2:3]   (y_pred[3,1]*I(3))[:, 3]    I(3)]


for i = 2:size(yVals)[2]
    row = [(y_pred[1,i]*I(3))   (y_pred[2,i]*I(3))[:, 2:3]   (y_pred[3,i]*I(3))[:, 3]    I(3)]

    global A = [A; row]
end


y_meas_vec = reshape(y_meas, length(y_meas));

params_guess = A \ y_meas_vec;

params = gauss_newton(params_guess)

T_hat = params_to_matrix(params[1:6]) # reshape(params[1:6], 3, 3)'
if T_hat[3,3] < 0
    # println("Correcting Sign Ambiguity in T₆!") # Should only be negative if the non-orthogonality angle is >90°, which should never happen
    T_hat[3,3] = -T_hat[3,3]
end

if T_hat[2,2] < 0
    # println("Correcting double sign ambiguits in T₄ and T₅!")
    T_hat[2,2] = -T_hat[2,2]
    T_hat[3,2] = -T_hat[3,2]
end

if T_hat[1,1] < 0
    println("Error! Scale factor cannot be negative!")
    # @error "Error! Scale factor cannot be negative!"
    scale_flag = true
end

bias_hat = params[7:end]

a_hat, b_hat, c_hat, ρ_hat, ϕ_hat, λ_hat = extractParameters(T_hat)
bx_hat, by_hat, bz_hat = bias_hat
ϵ_a, ϵ_b, ϵ_c, ϵ_ρ, ϵ_ϕ, ϵ_λ, ϵ_bx, ϵ_by, ϵ_bz = 
            (s_a - a_hat), (s_b - b_hat), (s_c - c_hat), (ρ - ρ_hat), (ϕ - ϕ_hat), (λ - λ_hat), (bias_x - bx_hat), (bias_y - by_hat), (bias_z - bz_hat)

err = sum(abs.(T_matrix - T_hat))

if err > 0.3 
    error_flag = true 
end


if error_flag || scale_flag
    println("\n\n")
    println("_________________________LINEAR TEST_________________________")
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
    println("\n\n")

    y_cor = (T_hat^(-1))*((y_meas) .- bias_hat)

    println(reshape(true_params[1:9], 3, 3)')
    println(T_hat)

else
    # println("Done")
    ttttt = 6
end



# plt1 = plot( y_meas[1,:], title = "Mag Fields (Before Correction)", xlabel = "Time", ylabel = "μT", color = :red, label = "Measured")
# plt1 = plot!(yᴮ_true[1,:], color = :red, linestyle = :dash, label = "Truth (?)")

# plt2 = plot(y_meas[2,:], color = :blue, label = false)
# plt2 = plot!(yᴮ_true[2,:], color = :blue, linestyle = :dash, label = false)

# plt3 = plot(y_meas[3,:], color = :green, label = false)
# plt3 = plot!(yᴮ_true[3,:], color = :green, linestyle = :dash, label = false)

# display(plot(plt1, plt2, plt3, layout = (3, 1)))



# pltAft1 = plot(yᴮ_true[1,:], title = "Mag Fields (After Correction)", color = :red, linestyle = :dash, label = "Truth (?)")
# pltAft1 = plot!(y_cor[1,:], xlabel = "Time", ylabel = "μT", color = :red, label = "Adjusted")

# pltAft2 = plot(yᴮ_true[2,:], color = :blue, linestyle = :dash, label = false)
# pltAft2 = plot!(y_cor[2,:], color = :blue, label = false)

# pltAft3 = plot(yᴮ_true[3,:], color = :green, linestyle = :dash, label = false)
# pltAft3 = plot!(y_cor[3,:], color = :green, label = false)

# display(plot(pltAft1, pltAft2, pltAft3, layout = (3, 1)))


# y_dif = abs.(y_meas - yᴮ_true)
# pltDif1 = plot( y_dif[1,:], color = :red, title = "Init Error")
# pltDif1 = plot!(y_dif[2,:], color = :blue, label = false)
# pltDif1 = plot!(y_dif[3,:], color = :green, label = false)

# y_dif = abs.(y_cor - yᴮ_true)
# pltDif2 = plot( y_dif[1,:], color = :red, title = "Final Error (60s)")
# pltDif2 = plot!(y_dif[2,:], color = :blue, label = false)
# pltDif2 = plot!(y_dif[3,:], color = :green, label = false)

# display(plot(pltDif1, pltDif2, layout = (2, 1)))