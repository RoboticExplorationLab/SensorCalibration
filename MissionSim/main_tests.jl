#########
# SETUP #
#########
using Infiltrator, Test

using StaticArrays, SatelliteDynamics, EarthAlbedo 
using Distributions, LinearAlgebra, Plots, JLD2, Random
using ProgressMeter


include("mag_field.jl");     # Contains IGRF13 stuff 
include("rotationFunctions.jl");

include("CustomStructs.jl");         using .CustomStructs
include("Simulator/Simulator.jl");   using .Simulator
include("Estimator/Estimator.jl");   using .Estimator 
include("Controller/Controller.jl"); using .Controller

@enum(Operation_mode, mag_cal = 1, detumble, diode_cal, mekf, chill, finished)
Base.to_index(om::Operation_mode) = Int(s)
include("state_machine.jl")


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
    x = [r₀; v₀; q₀; ω₀; β₀]

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

    refl = load_refl("data/refl.jld2", scale)  
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



""" Test 1 - MEKF 

  Run MEKF on ideal/noiseless data

  Expected Results:
    - ||Pchol|| < ||Pchol_pred||
    - ||Pchol|| > ||Pchol⁺||
    - ||x⁺ - xₜ|| < ||x₀ - xₜ||  (output is closer to truth than input)
"""
@testset "1) MEKF" begin 

    _Re = 6378136.3                 # Radius of the earth (m)) 
    _mu = (3.9860044188)*(10^(14))  # Standard Gravitational Parameter (m^3/s^2) 
    _E_am0 = 1366.9 # Irradiance of sunlight (TSI - visible & infrared), W/m^2 
    CONSTANTS = (_Re = _Re, _mu = _mu, _E_am0 = _E_am0)

    @testset "1a) Prediction" begin
        """ Test 1a - Prediction 
        Expected Results:
            - ||Pchol|| < ||Pchol_pred||
            - ||Pchol|| > ||Pchol⁺||
            - ||x⁺ - xₜ|| < ||x₀ - xₜ||  (output is closer to truth than input)

            (It is if we rotate ω before prediction...)
        """

        # Generate an orbit, and set our satellite to the true values 
        x₀, _ = get_initial_state(; detumbled = false, bias_less = false);
        # x₀ = STATE(x₀.r, x₀.v, x₀.q, SVector{3, Float64}(deg2rad(7), -deg2rad(10), deg2rad(10)), x₀.β)
        # x₀ = STATE(x₀.r, x₀.v, SVector{4, Float64}(1, 0, 0, 0), x₀.ω, x₀.β)

        mag = MAGNETOMETER(ones(3), zeros(3), zeros(3))
        dio = DIODES(ones(6), [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi], [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)])
        sat = SATELLITE( [0.2 0 0; 0 0.2 0; 0 0 0.2], mag, dio,  [x₀[7:10]; x₀[14:16]],  ones(6, 6))
        ω̃   = x₀[11:13] + x₀[14:16] # prediction subtracts this out, so we need to add it in 
        dt  = 1.0;
        N   = 6
        t   = Epoch(2020, 1, 1)
        u   = [0, 0, 0]

        # ### Plot SAT_STATE_PREDICTED vs actual STATE (Should be approximately equal for small dt)
        # dt = 0.2
        # Nn = 37
        # v₀ = [1, -1, 2]; v₀ /= norm(v₀)
        # sat_state = SAT_STATE(; q = x₀.q, β = x₀.β);
        # vᴮ = zeros(Nn, 3)
        # v̂ᴮ = zeros(Nn, 3)

        # xqs = zeros(Nn, 4)
        # sqs = zeros(Nn, 4)
        # x = x₀
        # for i = 1:Nn 
        #     # Predict 
        #     global sat_state, _ = Estimator.prediction(sat_state, ω̃, dt, N)
        #     ᴮQᴵ = (H' * L(sat_state.q) * R(sat_state.q)' * H)'
        #     v̂ᴮ[i, :]  .= ᴮQᴵ * v₀
        #     sqs[i, :] .= sat_state.q

        #     # Update        
        #     global x = rk4(sat.J, x, u, t, dt)
        #     vᴮ[i, :]  .= (H' * L(x.q) * R(x.q)' * H)' * v₀
        #     xqs[i, :] .= x.q

            
        #     # RESET so as to not accumulate error 
        #     # sat_state = SAT_STATE(; q = x.q, β = x.β);
        # end
        # plot(vᴮ, c = [:red :blue :green]); plot!(v̂ᴮ, c = [:red :blue :green], ls = :dash)
        # #########

        # Generate actual next state
        x⁺ = rk4(sat, x₀, u, t, dt)

        # Predict next state
        ᴵQᴮ = I(3) #dcm_from_q(x₀[7:10])    # (H' * L(sat_state.q) * R(sat_state.q)' * H)
        sat_state_pred, A = Estimator.prediction(sat.state, ᴵQᴮ * ω̃ , dt)

        # Verify that our prediction is closer than our initial guess 
        function qErr(q₁, q₂)
            return norm( qmult( qconj(q₁), q₂) - [0, 0, 0, 1])
        end
        qErr₀ = qErr(sat.state[1:4],  x⁺[7:10]) # norm(cayley_map(sat.state.q, x⁺.q))
        qErr₁ = qErr(sat_state_pred[1:4], x⁺[7:10]) # norm(cayley_map(sat_state_pred.q, x⁺.q))
        @test qErr₀ ≥ qErr₁
    end 

    @testset "1b) Mag Measurement" begin 
        """ Test 1b - Mag Measurement 
        Process:
          Pass in true x and noiseless mag measurement 

        Expected Results:
            - Bᴮ_exp == B̃ᴮ
        """

        # Generate an orbit, and set our satellite to the true values 
        x₀, _ = get_initial_state(; detumbled = false, bias_less = false);
        mag = MAGNETOMETER(ones(3), zeros(3), zeros(3))
        dio = DIODES(ones(6), [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi], [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)])
        sat = SATELLITE( [0.2 0 0; 0 0.2 0; 0 0 0.2], mag, dio,  [x₀[7:10]; x₀[14:16]],  ones(6, 6))
       
        dt  = 1.0;
        N   = 6
        t   = Epoch(2020, 1, 1)
        ᴮQᴵ = dcm_from_q(x₀[7:10])'

        Bᴵ, Bᴮ, B̃ᴮ = Simulator.generate_magnetic_field(x₀[1:3], t, sat, ᴮQᴵ, dt; σ = 0.0)
        Bᴮ_exp, Hb = Estimator.mag_measurement(sat.state, Bᴵ)
        @test Bᴮ ≈ Bᴮ_exp 
        @test correct_mag_field(sat, B̃ᴮ) ≈ Bᴮ_exp 
    end

    @testset "1c) Diode Measurement" begin 
        """ Test 1b - Mag Measurement 
        Process:
          Pass in true x and noiseless mag measurement 

        Expected Results:
            - Bᴮ_exp == B̃ᴮ
        """

        # Generate an orbit, and set our satellite to the true values 
        x₀, _ = get_initial_state(; detumbled = false, bias_less = false);
        mag = MAGNETOMETER(ones(3), zeros(3), zeros(3))
        dio = DIODES(ones(6), [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi], [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)])
        sat = SATELLITE( [0.2 0 0; 0 0.2 0; 0 0 0.2], mag, dio,  [x₀[7:10]; x₀[14:16]],  ones(6, 6))

        dt  = 1.0;
        N   = 6
        t   = Epoch(2020, 1, 1)
        alb = get_albedo(2)

        _Re = 6378136.3                 # Radius of the earth (m)) 
        _mu = (3.9860044188)*(10^(14))  # Standard Gravitational Parameter (m^3/s^2) 
        _E_am0 = 1366.9 # Irradiance of sunlight (TSI - visible & infrared), W/m^2 
        CONSTANTS = (_Re = _Re, _mu = _mu, _E_am0 = _E_am0)

        ᴮQᴵ = dcm_from_q(x₀[7:10])'
        sᴵ, sᴮ, ecl = Simulator.update_sun_vectors(x₀[1:3], t, ᴮQᴵ, dt)
        It, Ĩ, ηI  = Simulator.generate_diode_currents(sat, x₀[1:3], alb, sᴵ, sᴮ, ecl, CONSTANTS; σ = 0.0)
        C, α, ϵ    = sat.diodes.calib_values, sat.diodes.azi_angles, sat.diodes.elev_angles
        I_exp,  _  = Estimator.current_measurement(sat.state, C, α, ϵ, sᴵ, 6,  x₀[1:3], t, alb)

        @test It ≈ I_exp 
        @test Ĩ  ≈ I_exp

    end

    # DOESN"T work with initial bias...? Got mag meas working...
    @testset "1 Together" begin 
        """
             Verify that our predicted state matches the true next state better than 
            our original guess did. Also verify it leads to better mag and diode measurements

            (All noise less)
        """
        # Starting point 
        # Rk4 + Update sensor measurements (no noise) 
        # Predict next x 
        # Use that next x to predict mag measurements (should be better than initial guess)

        x₀, _ = get_initial_state(; detumbled = false, bias_less = false);
        mag = MAGNETOMETER(ones(3), zeros(3), zeros(3))
        dio = DIODES(ones(6), [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi], [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)])
        sat = SATELLITE( [0.2 0 0; 0 0.2 0; 0 0 0.2], mag, dio,  [x₀[7:10]; x₀[14:16]],  ones(6, 6))

        function qErr(q₁, q₂)
            return norm( qmult( qconj(q₁), q₂) - [0, 0, 0, 1])
        end

        dt  = 0.2;  # Needs to be somewhat small
        N   = 6
        t   = Epoch(2020, 1, 1)
        u   = SVector{3, Float64}(0, 0, 0)
        alb = get_albedo(2)

        # Update state and measurements
        x⁺ = rk4(sat, x₀, u, t, dt)
        sat⁺ = SATELLITE(sat.J, sat.magnetometer, sat.diodes, [x⁺[7:10]; x⁺[14:16]], sat.covariance)

        truth, sens, ecl, noise = Simulator.generate_measurements(Simulator.SIM(0.0), sat⁺, alb, x⁺, t, CONSTANTS, dt) #; σB = 0.0, σ_gyro_scale = 0.0, σr = 0.0, σ_current_scale = 0.0)
        sensors = correct_mag_field(sat⁺, sens.magnetometer)
        @debug sens.magnetometer == sensors # If no bias 

        # Predict next state
        # GYRO neeeds to be rotated to inertial...? But NOT β...
        ᴵQᴮ =  I(3) #dcm_from_q(x₀[7:10])   # SMatrix{3, 3, Float64, 9}(I(3)) #(H' * L(sat_state₀.q) * R(sat_state₀.q)' * H) # TRUE rotation
        sat_state_pred, A = Estimator.prediction(sat.state, ᴵQᴮ * (sens.gyro - x₀[14:16]) + x₀[14:16], dt)

        # A[1:3, 4:6] .= ᴵQᴮ' * A[1:3, 4:6]  # Unrotate the bias portion back to body frame (IDK...)
        # A[4:6, 4:6] .= ᴵQᴮ' * A[4:6, 4:6]

        #   Verify prediction is closer to truth than initial guess 
        qErr₀ = qErr(sat.state[1:4],      x⁺[7:10]) # norm(cayley_map(sat.state.q, x⁺.q))
        qErr₁ = qErr(sat_state_pred[1:4], x⁺[7:10]) # norm(cayley_map(sat_state_pred.q, x⁺.q))
        @test qErr₀ ≥ qErr₁

        # Use predicted state to predict mag measurements 
        Bᴮ_exp⁺, Hb⁺ = Estimator.mag_measurement(sat_state_pred, truth.Bᴵ_hist)
        Bᴮ_exp₀, Hb₀ = Estimator.mag_measurement(sat.state,  truth.Bᴵ_hist)
        @test norm(Bᴮ_exp⁺ - truth.Bᴮ_hist) < norm(Bᴮ_exp₀ - truth.Bᴮ_hist) && (norm(Bᴮ_exp⁺ - truth.Bᴮ_hist) > 0)

        # Use predicted state to predict diode measurements 
        C, α, ϵ    = sat.diodes.calib_values, sat.diodes.azi_angles, sat.diodes.elev_angles
        sᴵ = truth.sᴵ_hist / norm(truth.sᴵ_hist)
        I_exp⁺,  Hc⁺  = Estimator.current_measurement(sat_state_pred, C, α, ϵ, sᴵ, 6, x₀[1:3], t, alb)
        I_exp₀,  Hc₀  = Estimator.current_measurement(sat.state,      C, α, ϵ, sᴵ, 6, x₀[1:3], t, alb)
        # @test norm(I_exp⁺ - truth.I_hist) < norm(I_exp₀ - truth.I_hist) && (norm(I_exp⁺ - truth.I_hist) > 0)
        @test norm(I_exp⁺ - sens.diodes) < norm(I_exp₀ - sens.diodes) && (norm(I_exp⁺ - sens.diodes) > 0)
    

        # Innovation 
        z⁺ = [sens.magnetometer - Bᴮ_exp⁺; sens.diodes - I_exp⁺];
        z₀ = [sens.magnetometer - Bᴮ_exp₀; sens.diodes - I_exp₀];
        H⁺ = [Hb⁺; Hc⁺];
        H₀ = [Hb₀; Hc₀];
        @test (norm(z⁺) < norm(z₀)) && (norm(z⁺) > 0)


        # Verify that our MEKF-updated guess is closer to the truth than our original point
        data = Estimator.MEKF_DATA();
        Pchol_yy = Estimator.qrᵣ([A * H⁺'; data.Vchol]);    
        Ll = (((A' * A * H⁺') / Pchol_yy) / Pchol_yy');
        @test norm(Ll * z⁺) < norm(Ll * z₀)# Not actually strictly necessary 
        state⁺, diodes⁺ = Estimator.update(sat_state_pred, sat.diodes, Ll, z⁺; calibrate_diodes = false)
        # state₀, diodes₀ = Estimator.update(sat_state_pred, sat.diodes, Ll, z₀; calibrate_diodes = false)
        eq⁺ = norm( L(state⁺.q)' * x⁺.q - [1, 0, 0, 0])
        eq₀ = norm( L(sat_state₀.q)' * x⁺.q - [1, 0, 0, 0])
        eβ⁺ = norm( state⁺.β - x⁺.β)
        eβ₀ = norm( sat_state₀.β - x⁺.β)
        @test eq⁺ < eq₀
        @test eβ⁺ < eβ₀

        # Do it all together now... (with no noise)
        Pchol = sat.covariance[1:6, 1:6]
        x̂⁺, d̂⁺, P̂⁺ = Estimator.sqrt_mekf(sat.state, sat.diodes, Pchol, alb, sensors.gyro, truth.Bᴵ, truth.sᴵ, 
                                          sensors.magnetometer, sensors.diodes, sensors.pos, dt, UpperTriangular(data.Wchol[1:6, 1:6]), data.Vchol; calibrate_diodes = false)
        
        
    end
end


#######################
# TEST FRAME OF OMEGA #
#######################
@testset "Omega frame" begin
    x₀ = SVector{16, Float64}([
        [3e6; 3.8e6; -4.8e6];   # r 
        [-6861; 2331; -2476];   # v 
        [[0, 1, 0] * sin(deg2rad(45)); cos(deg2rad(45))];   # q 
        [deg2rad(10.0); deg2rad(0.0); deg2rad(0.0)];        # ω
        zeros(3)
    ])
    # x₀, _ = get_initial_state(; detumbled = true, bias_less = true);
    # # q = randn(4); q /= norm(q)
    # x₀ = SVector{16, Float64}([x₀[1:6];   # r, v 
    #                           [cos(deg2rad(45)); [0, 1, 0] * sin(deg2rad(45))];   # q
    #                           [deg2rad(0.0), deg2rad(0.0), deg2rad(10.0)];       # ω
    #                           x₀[14:16]]  )                                      # β 

    mag = MAGNETOMETER(ones(3), zeros(3), zeros(3))
    dio = DIODES(ones(6), [0.0;      pi;     (pi/2); (-pi/2); 0.0;    pi], [(-pi/4); (pi/4);  0.0;    0.0;    (pi/4); (-pi/4)])
    sat = SATELLITE( [0.2 0 0; 0 0.2 0; 0 0 0.2], mag, dio,  [x₀[7:10]; x₀[14:16]],  ones(6, 6))

    vᴵ = [1, 0.0, 0.0]
    J  = [1.0 0 0; 0 1 0; 0 0 1]
    t  = Epoch(2020, 1, 1)
    u  = [0.0, 0.0, 0.0]
    dt = 0.1

    N = 360 
    vᴮ = zeros(N, 3);
    v̂ᴮ = zeros(N, 3);
    qs = zeros(N, 4);
    q̂s = zeros(N, 4);
    ωs = zeros(N, 3);
    x = x₀;
    ss = [x₀[7:10]; x₀[14:16]];
    for i = 1:N
        ᴵQᴮ = dcm_from_q(x[7:10])
        vᴮ[i, :] .= ᴵQᴮ' * vᴵ;
        qs[i, :] .= x[7:10]
        global ss, _ = Estimator.prediction(ss, x[11:13], dt)
        ᴵQ̂ᴮ = dcm_from_q(ss[1:4])
        v̂ᴮ[i, :] .= ᴵQ̂ᴮ' * vᴵ;
        q̂s[i, :] .= ss[1:4]
        ωs[i, :] .= x[11:13]
        global x = rk4(sat, x, u, t, dt);
    end

    plot( [0, 1], [0, 0], [0, 0], label = "+X", c = :red,   lw = 3);
    plot!([0, 0], [0, 1], [0, 0], label = "+Y", c = :blue,  lw = 3);  
    plot!([0, 0], [0, 0], [0, 1], label = "+Z", c = :green, lw = 3, 
                xlim = [-1.1, 1.1], ylim = [-1.1, 1.1], zlim = [-1.1, 1.1] );
    @gif for i ∈ 1:N
        plot!([0, vᴮ[i, 1]], [0, vᴮ[i, 2]], [0, vᴮ[i, 3]], c = :violet, lw = 2, label = false)
                    # xlim = [-1.1, 1.1], ylim = [-1.1, 1.1], zlim = [-1.1, 1.1] )
    end every 5

    plot(vᴮ, title = "RK4", label = ["x" "y" "z"], c = [:red :blue :green])
    plot(v̂ᴮ, title = "Pred", label = ["x" "y" "z"], c = [:red :blue :green])

end




#########################
# Tests for Calibrating #
#########################
# -> Gyro and bias match true, but we are trying to get C, α, and ϵ
# There will be no position, magnetometer, or gyro noise/bias, but there will be current noise (and the angle maybe larger...)

@testset "Test for Calibrating" begin 
    x₀, _ = get_initial_state(; detumbled = true, bias_less = true)
    # x₀ = STATE(x₀.r, x₀.v, SVector{4, Float64}([cos(deg2rad(45)); [0, 1, 0] * sin(deg2rad(45))] ), SVector{3, Float64}(deg2rad(10.0), deg2rad(0.0), deg2rad(0.0)), x₀.β)
    J  = SMatrix{3, 3, Float64, 9}([1.0 0 0; 0 1.0 0; 0 0 1.0])
    t  = Epoch(2020, 1, 1)
    u  = SVector{3, Float64}(0.0, 0.0, 0.0)
    dt = 0.5

    sat_state₀ = SAT_STATE(x₀.q, x₀.β);
    sat_truth = SATELLITE(; J = J, mag = MAGNETOMETER(; ideal = true), dio = DIODES(), sta = sat_state₀);
    sat_est   = SATELLITE(; J = J, mag = MAGNETOMETER(; ideal = true), dio = DIODES(; ideal = true), sta = sat_state₀);
    alb       = get_albedo(2) ;
    noise     = MEKF_DATA();

    # W = Diagonal( [0.001 * ones(3); 0.001 * ones(3); 0.1 * ones(6); deg2rad(3) * ones(6); deg2rad(3) * ones(6)].^2 )
    # Wchol = Estimator.chol(Matrix(W)) 
    # V = Diagonal( [0.001 * ones(3); 0.001 * ones(6)].^2 )
    # Vchol = Estimator.chol(Matrix(V)) 
    # noise     = MEKF_DATA(Wchol, Vchol)

    N = 3000
    x = x₀;
    Cs = zeros(N, 6); αs = zeros(N, 6); ϵs = zeros(N, 6);
    ecls = zeros(N);
    # rs = zeros(N, 3); vs = zeros(N, 3); qs = zeros(N, 4); ωs = zeros(N, 3); βs = zeros(N, 3);  

    for i = 1:N
        Cs[i, :] .= sat_est.diodes.calib_values
        αs[i, :] .= rad2deg.(sat_est.diodes.azi_angles)
        ϵs[i, :] .= rad2deg.(sat_est.diodes.elev_angles)
        # rs[i, :] .= x.r 
        # vs[i, :] .= x.v 
        # qs[i, :] .= x.q 
        # ωs[i, :] .= x.ω 
        # βs[i, :] .= x.β
        global x = rk4(J, x, u, t, dt; σβ = 0.0);
        global sat_truth  = SATELLITE(sat_truth.J, sat_truth.magnetometer, sat_truth.diodes, SAT_STATE(x.q, x.β), sat_truth.covariance)
        tr, sensors, ecl, nse = generate_measurements(sat_truth, alb, x, t, dt; σB = deg2rad(0.0), σr = 0.0, σ_gyro_scale = 0.0, σ_current_scale = 0.01)
        ecls[i] = ecl
        if ecl > 0.1
            global sat_est = Estimator.estimate(sat_est, sensors, noise, alb, t, dt; calibrate_diodes = true)
            st = SAT_STATE(x.q, x.β)
            global sat_est = SATELLITE(sat_est.J, sat_est.magnetometer, sat_est.diodes, st, sat_est.covariance)
        end
        global t = t + dt 
    end;

    cPlt = []; αPlt = []; ϵPlt = []
    for i = 1:6 
        cp = plot(   Cs[:, i], label = false, c = :red);
        cp = hline!([sat_truth.diodes.calib_values[i]], ls = :dash, c = :green, label = false); #"Truth");
        cp = hline!([1.0], ls = :dot, c = :blue, label = false);# "Init Guess");
        cp = plot!(ylim = [0.8, 1.2])
        push!(cPlt, cp);

        ap = plot(   αs[:, i], label = false, c = :red);
        ap = hline!([rad2deg(sat_truth.diodes.azi_angles[i])], ls = :dash, c = :green, label = false);
        ap = hline!([αs[1, i]], ls = :dot, c = :blue, label = false);
        ap = plot!( ylim = [rad2deg(sat_truth.diodes.azi_angles[i]) - 5, rad2deg(sat_truth.diodes.azi_angles[i]) + 5])
        push!(αPlt, ap);

        ep = plot(   ϵs[:, i], label = false, c = :red);
        ep = hline!([rad2deg(sat_truth.diodes.elev_angles[i])], ls = :dash, c = :green, label = false);
        ep = hline!([ϵs[1, i]], ls = :dot, c = :blue, label = false);
        ep = plot!( ylim = [rad2deg(sat_truth.diodes.elev_angles[i]) - 5, rad2deg(sat_truth.diodes.elev_angles[i]) + 5] )
        push!(ϵPlt, ep);
    end;
    display(plot(cPlt..., title = "C", layout = [2, 3]))
    display(plot(αPlt..., title = "α", layout = [2, 3]))
    display(plot(ϵPlt..., title = "ϵ", layout = [2, 3]))

end;









# @testset "Add to dynamics tests" begin 
#     # Does my RK4 rotate in a circle correctly?
#     q = randn(4); q /= norm(q)
#     q₀ = SVector{4, Float64}(q) #(1, 0, 0, 0)
#     β  = SVector{3, Float64}(0, 0, 0)
#     ω  = SVector{3, Float64}(0.0, deg2rad(10), 0.0)
#     J  = SMatrix{3, 3, Float64, 9}([0.2 0 0; 0 0.2 0; 0 0 0.2])
#     u  = SVector{3, Float64}(0.0, 0, 0)

#     x = STATE(; q = q₀, β = β, ω = ω)
#     t = Epoch(2020, 1, 1)
#     dt = 1.0
#     N = Int(36 / dt)
#     qs = zeros(N, 4)
#     for i = 1:N
#         qs[i, :] .= x.q
#         global x = Simulator.rk4(J, x, u, t, dt; σβ = 0.0)  # No bias
#     end
#     # display(plot(qs))

#     @test norm(cayley_map(x.q, q₀)) ≈ 0.0 atol = 1e-5
# end



function rotDyn(J, x)
    q = x[1:4]
    ω = x[5:7]

    Lq = [q[1] -q[2:4]';  q[2:4] (q[1]*I(3)+hat(q[2:4])) ]
    q̇ = 0.5 * Lq * [0; ω]

    ω̇ = J \ (-cross(ω, J * ω))

    return [q̇; ω̇ ]
end

function rk4(J, x, dt)
    k₁ = dt * rotDyn(J, x)
    k₂ = dt * rotDyn(J, x + k₁/2)
    k₃ = dt * rotDyn(J, x + k₂/2)
    k₄ = dt * rotDyn(J, x + k₃)

    x⁺ = x + (1/6) * (k₁ + 2 * k₂ + 2 * k₃ + k₄)
    x⁺[1:4] /= norm(x⁺[1:4])
    return x⁺
end


# Check and visualize ω
q₀ = randn(4); q₀ /= norm(q₀)
ω₀ = [deg2rad(10.0), deg2rad(0.0), deg2rad(0.0)]
x₀ = [q₀; ω₀]
J  = [0.2 0 0; 0 0.2 0; 0 0 0.2]
dt = 1.0 
N  = 36 

v₀ = [0, 0, 1]

x = x₀
vᴮ = zeros(3, N + 1)
for i = 1:N + 1
    ᴵQᴮ =  H' * L(x[1:4]) * R(x[1:4])' * H 
    @show ᴵQᴮ'
    vᴮ[:, i] .= ᴵQᴮ' * v₀
    global x = rk4(J, x, dt)
end

plot( [0, 1], [0, 0], [0, 0], label = "+X", c = :red, lw = 3);
plot!([0, 0], [0, 1], [0, 0], label = "+Y", c = :blue, lw = 3);  
plot!([0, 0], [0, 0], [0, 1], label = "+Z", c = :green, lw = 3, 
            xlim = [-1.1, 1.1], ylim = [-1.1, 1.1], zlim = [-1.1, 1.1] );
@gif for i ∈ 1:N + 1
    plot!([0, vᴮ[1, i]], [0, vᴮ[2, i]], [0, vᴮ[3, i]], c = :violet, lw = 2, label = false)
                # xlim = [-1.1, 1.1], ylim = [-1.1, 1.1], zlim = [-1.1, 1.1] )
end every 1

# scatter(vᴮ[3, :], vᴮ[1, :])
# plot(vᴮ')








# Check and visualize q update 
q₀ = [0, 1, 1, 0]; q₀ /= norm(q₀)
ω₀ = 0.1 * randn(3) #[deg2rad(10.0), deg2rad(-2.0), deg2rad(3.0)]
x₀ = [q₀; ω₀]
J  = [0.2 0 0; 0 0.2 0; 0 0 0.2]
dt = 0.2
N  = Int(36 / dt) 

x = x₀
q̂s = zeros(N + 1, 4)
qs = zeros(N + 1, 4)
for i = 1:N + 1
    ᴵQᴮ = (H' * L(x[1:4]) * R(x[1:4])' * H) # TRUE rotation
        
    ω = ᴵQᴮ * x[5:7];
    r = ω / norm(ω);
    θ = dt * norm(ω); 
    q⁺ = qmult(x[1:4], [cos(θ / 2); r*sin(θ / 2)]);

    q̂s[i, :] .= q⁺; 


    global x = rk4(J, x, dt);
    qs[i, :] .= x[1:4];

end

plot(qs, c = [:red :blue :green :violet], ls = :dash);
plot!(q̂s, c = [:red :blue :green :violet])



