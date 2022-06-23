#########
# SETUP #
#########


path = "../src/MissionSim/"
##### Stuff needed to run in isolation ##########################################
# using Infiltrator, Test
# using StaticArrays, SatelliteDynamics, EarthAlbedo 
# using Distributions, LinearAlgebra, Plots, JLD2, Random
# using ProgressMeter

# include(path * "mag_field.jl");     # Contains IGRF13 stuff 
# include(path * "quaternions.jl")

# include(path * "CustomStructs.jl");         using .CustomStructs
# include(path * "Simulator/Simulator.jl");   using .Simulator
# include(path * "Estimator/Estimator.jl");   using .Estimator 
# include(path * "Controller/Controller.jl"); using .Controller

# @enum(Operation_mode, mag_cal = 1, detumble, diode_cal, mekf, chill, finished)
# Base.to_index(om::Operation_mode) = Int(s)
# include(path * "state_machine.jl")
#################################################################################


function get_initial_state(; _Re = 6378136.3, detumbled = false, bias_less = false) 
    inc = rand(1:180);  # 51.6426 
    Ω   = rand(1:360);  # 178.1369 
    ω   = rand(1:360);  # 174.7410 
    M   = rand(1:360);  # 330.7918  (+94/95 is just before sun, -40 is just before eclipse)
    z   = rand(300:1000) * 1e3;  # 421e3 (Altitude)

    ecc_max = (z / _Re)  # SemiMajor Axis = Apogee / (1 + ecc), which needs to be > _Re
    ecc = rand(0:0.001:0.025); #ecc_max);       # 0.0001717 

    sma = (_Re + z) / (1 + ecc);  # Apogee = semi_major * (1 + ecc)

    oe0 = [sma, ecc, inc, Ω, ω, M];   # Initial state, oscullating elements
    eci0 = sOSCtoCART(oe0, use_degrees = true); # Convert to Cartesean

    r₀ = SVector{3, Float64}(eci0[1:3]);
    v₀ = SVector{3, Float64}(eci0[4:6]);
    q₀ = randn(4);  q₀ = SVector{4, Float64}(q₀ / norm(q₀));
    ω₀ = (detumbled) ? SVector{3, Float64}(0.05 * randn(3)) : SVector{3, Float64}(0.5 * randn(3));
    β₀ = (bias_less) ? SVector{3, Float64}(0.0, 0.0, 0.0)  : SVector{3, Float64}(rand(Normal(0.0, deg2rad(2.0)), 3));
    
    T_orbit = orbit_period(oe0[1]);
    x = STATE(r₀, v₀, q₀, ω₀, β₀);

    return x, T_orbit 
end;

function get_albedo(scale = 1) 

    function load_refl(path = "data/refl.jld2", scale = 1)
        temp = load(path)
    
        refl = REFL( temp["data"][1:scale:end, 1:scale:end], temp["type"], temp["start_time"], temp["stop_time"])
    
        return refl
    end
    lat_step = 1.0 * scale
    lon_step = 1.25 * scale

    refl = load_refl(path * "data/refl.jld2", scale)  
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



# """ Test 1 - MEKF 

#   Run MEKF on ideal/noiseless data

#   Expected Results:
#     - ||Pchol|| < ||Pchol_pred||
#     - ||Pchol|| > ||Pchol⁺||
#     - ||x⁺ - xₜ|| < ||x₀ - xₜ||  (output is closer to truth than input)
# """
@testset "1) MEKF" begin 

    function qErr(q₁, q₂)
        return norm((L(q₁)' * q₂)[2:4])  # Equivalent to subtracting [±1, 0, 0, 0], but easier
    end


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
        sat_state = SAT_STATE(; q = x₀.q, β = x₀.β);
        sat = SATELLITE(; J = SMatrix{3, 3, Float64, 9}([0.2 0 0; 0 0.2 0; 0 0 0.2]), sta = sat_state);
        ω̃   = x₀.ω + x₀.β # prediction subtracts this out, so we need to add it in 
        dt  = 1.0;
        N   = 6
        t   = Epoch(2020, 1, 1)
        u   = SVector{3, Float64}(0, 0, 0)

        # Generate actual next state
        x⁺ = rk4(sat.J, x₀, u, t, dt)

        # Predict next state
        sat_state_pred, A = Estimator.prediction(sat.state, ω̃ , dt, N)

        # Verify that our prediction is closer than our initial guess 

        qErr₀ = qErr(sat.state.q,      x⁺.q) # norm(cayley_map(sat.state.q, x⁺.q))
        qErr₁ = qErr(sat_state_pred.q, x⁺.q) # norm(cayley_map(sat_state_pred.q, x⁺.q))
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
        sat_state = SAT_STATE(; q = x₀.q, β = x₀.β);
        mag = MAGNETOMETER(ones(3), zeros(3), x₀.β)
        sat = SATELLITE(; J = SMatrix{3, 3, Float64, 9}([0.2 0 0; 0 0.2 0; 0 0 0.2]), sta = sat_state, mag = mag, dio = DIODES(; ideal = true));
        
        # sat = SATELLITE(; J = SMatrix{3, 3, Float64, 9}([0.2 0 0; 0 0.2 0; 0 0 0.2]), sta = sat_state, mag = MAGNETOMETER(; ideal = true));
        dt  = 1.0;
        N   = 6
        t   = Epoch(2020, 1, 1)
        ᴮQᴵ = quat2rot(x₀.q)'

        Bᴵ, Bᴮ, B̃ᴮ = Simulator.mag_measurement(sat, x₀, ᴮQᴵ, t, dt; σ = 0.0)
        Bᴮ_exp, Hb = Estimator.mag_measurement(sat.state, Bᴵ, N; calibrate_diodes = false)
        @test Bᴮ ≈ Bᴮ_exp 
        @test correct_magnetometer(sat, B̃ᴮ) ≈ Bᴮ_exp 
    end

    @testset "1c) Diode Measurement" begin 
        """ Test 1b - Mag Measurement 
        Process:
          Pass in true x and noiseless mag measurement 

        Expected Results:
            - Bᴮ_exp == B̃ᴮ
        """

        for i = 1:10
            # Generate an orbit, and set our satellite to the true values 
            x₀, _ = get_initial_state(; detumbled = false, bias_less = false);
            sat_state = SAT_STATE(; q = x₀.q, β = x₀.β);
            sat = SATELLITE(; J = SMatrix{3, 3, Float64, 9}([0.2 0 0; 0 0.2 0; 0 0 0.2]), sta = sat_state, mag = MAGNETOMETER(; ideal = true));
            dt  = 1.0;
            N   = 6
            t   = Epoch(2020, 1, 1)
            alb = get_albedo(2)

            ᴮQᴵ = quat2rot(x₀.q)'
            sᴵ, sᴮ, ecl = Simulator.sun_measurement( x₀, ᴮQᴵ, t)

            (ecl < 0.1) && continue 

            It, Ĩ, ηI   = Simulator.diode_measurement(sat, alb, x₀, ecl, sᴵ, sᴮ; σ_scale = 0.0, use_albedo = true)
            I_exp,  _   = Estimator.current_measurement(sat.state, sat.diodes, sᴵ, x₀.r, alb; use_albedo = true, calibrate_diodes = false)

            @test (It / norm(It)) ≈ I_exp 
            @test (Ĩ / norm(Ĩ))  ≈ I_exp
        end
    end

    # DOESN"T work with initial bias...? Got mag meas working...
    @testset "1 Together" begin 
        """
             Verify that our predicted state matches the true next state better than 
            our original guess did. Also verify it leads to better mag and diode measurements

            (All noise less)
        """

        for i = 1:100
            x₀, _ = get_initial_state(; detumbled = false, bias_less = true);
            sat_state₀ = SAT_STATE(; q = x₀.q, β = x₀.β);
            mag = MAGNETOMETER(ones(3), zeros(3), x₀.β)
            sat = SATELLITE(; J = SMatrix{3, 3, Float64, 9}([0.2 0 0; 0 0.2 0; 0 0 0.2]), sta = sat_state₀, mag = mag, dio = DIODES(; ideal = true));
            
            dt  = 0.2;  # Needs to be somewhat small
            N   = 6
            t   = Epoch(2020, 1, 1) + rand(1:365) * 24 * 60 * 60 
            u   = SVector{3, Float64}(0, 0, 0)
            alb = get_albedo(2);

            # Update state and measurements
            x⁺ = rk4(sat.J, x₀, u, t, dt)
            sat_state⁺ = SAT_STATE(; q = x⁺.q, β = x⁺.β);
            sat⁺ = SATELLITE(; J = sat.J, sta = sat_state⁺, dio = sat.diodes, mag = sat.magnetometer);
            truth, sensors, ecl, noise = generate_measurements(sat⁺, alb, x⁺, t, dt; σB = 0.0, σ_gyro_scale = 0.0, σr = 0.0, σ_current_scale = 0.0);

            data  = Estimator.MEKF_DATA()
            Pchol = sat.covariance[1:6, 1:6]

            # @debug sens == sensors # If no bias 
            @test x⁺.ω + x⁺.β == sensors.gyro

            # Predict next state
            sat_state_pred, A = Estimator.prediction(sat_state₀, sensors.gyro, dt, N; calibrate_diodes = false) 
            Pchol_pred = Estimator.qrᵣ([Pchol * A'; UpperTriangular(data.Wchol[1:6, 1:6])])

            #   Verify prediction is closer to truth than initial guess 
            qErr₀ = norm(cayley_map(sat.state.q, x⁺.q))
            qErr₁ = norm(cayley_map(sat_state_pred.q, x⁺.q))
            @test qErr₀ ≥ qErr₁



            # Use predicted state to predict mag measurements 
            Bᴮ_exp₀, Hb₀ = Estimator.mag_measurement(sat_state₀,     normalize(truth.Bᴵ), N; calibrate_diodes = false)
            Bᴮ_exp⁺, Hb⁺ = Estimator.mag_measurement(sat_state_pred, normalize(truth.Bᴵ), N; calibrate_diodes = false)
            @test norm(Bᴮ_exp⁺ - truth.Bᴮ) < norm(Bᴮ_exp₀ - truth.Bᴮ) && (norm(Bᴮ_exp⁺ - truth.Bᴮ) > 0)

            # Use predicted state to predict diode measurements 
            if ecl > 0.1
                I_exp₀,  Hc₀  = Estimator.current_measurement(sat_state₀,     sat.diodes, truth.sᴵ, sensors.pos, alb; use_albedo = true, calibrate_diodes = false)
                I_exp⁺,  Hc⁺  = Estimator.current_measurement(sat_state_pred, sat.diodes, truth.sᴵ, sensors.pos, alb; use_albedo = true, calibrate_diodes = false)
                @test norm(I_exp⁺ - normalize(truth.I)) < norm(I_exp₀ - normalize(truth.I)) && (norm(I_exp⁺ - normalize(truth.I)) > 0)

                # Innovation 
                z⁺ = [normalize(sensors.magnetometer) - Bᴮ_exp⁺; normalize(sensors.diodes) - I_exp⁺];
                
                ### This may get fixed, but for now...
                Hb⁺ *= 2 
                Hc⁺[:, 1:3] *= 2
                
                H⁺ = [Hb⁺; Hc⁺];


                # Verify that our MEKF-updated guess is closer to the truth than our original point

                Pchol_yy = Estimator.qrᵣ([Pchol_pred * H⁺'; data.Vchol]);

                
                L = (((Pchol_pred' * Pchol_pred * H⁺') / Pchol_yy) / Pchol_yy');
                state⁺, diodes⁺ = Estimator.update(sat_state_pred, sat.diodes, L, z⁺; calibrate_diodes = false)

                eq⁺ = qErr(state⁺.q, x⁺.q) 
                eq₀ = qErr(sat_state₀.q, x⁺.q)
                eβ⁺ = norm( state⁺.β - x⁺.β)
                eβ₀ = norm( sat_state₀.β - x⁺.β)
                @test eq⁺ < eq₀
                @test eβ⁺ < eβ₀


                # Do it all together now...
                x̂⁺, d̂⁺, P̂⁺ = Estimator.sqrt_mekf(sat.state, sat.diodes, Pchol, alb, sensors.gyro, truth.Bᴵ, truth.sᴵ, 
                                                sensors.magnetometer, sensors.diodes, sensors.pos, dt, UpperTriangular(data.Wchol[1:6, 1:6]), data.Vchol; calibrate_diodes = false)
                
                eq⁺ = qErr(x̂⁺.q, x⁺.q)
                eq₀ = qErr(sat_state₀.q, x⁺.q)
                eβ⁺ = norm( x̂⁺.β - x⁺.β)
                eβ₀ = norm( sat_state₀.β - x⁺.β)
                @test eq⁺ < eq₀
                @test eβ⁺ < eβ₀
            end
        end
    end
end


@testset "Diode Cal" begin
    @testset "Covariance" begin
        ### Setup 
        N   = 6
        dt  = 0.2;  # Needs to be somewhat small
        u   = SVector{3, Float64}(0, 0, 0)
        alb = get_albedo(2)
        
        data = MEKF_DATA() 

        Random.seed!()
        for i = 1:10 
            x₀, _ = get_initial_state(; detumbled = false, bias_less = true);
            sat_state₀ = SAT_STATE(; q = x₀.q, β = x₀.β);
            mag = MAGNETOMETER(ones(3), zeros(3), x₀.β)  ;
            sat = SATELLITE(; J = SMatrix{3, 3, Float64, 9}([0.2 0 0; 0 0.2 0; 0 0 0.2]), sta = sat_state₀, mag = mag,
                              dio = DIODES(; ideal = true));
            
            t   = Epoch(2020, 1, 1) + (rand(1:365) * 60 * 60 * 24);
            
            # Update state and measurements
            x⁺ = rk4(sat.J, x₀, u, t, dt);
            sat_state⁺ = SAT_STATE(; q = x⁺.q, β = x⁺.β);
            sat⁺ = SATELLITE(; J = sat.J, sta = sat_state⁺, dio = sat.diodes, mag = sat.magnetometer);

            truth, sens, ecl, noise = generate_measurements(sat⁺, alb, x⁺, t, dt; σB = 0.0, σ_gyro_scale = 0.0, σr = 0.0, σ_current_scale = 0.0);

            ecl
            # During an eclipse, covariance 𝑆ℎ𝑜𝑢𝑙𝑑 increase/do weird things
            if ecl < 0.3
                continue
            end

            sensors = correct_magnetometer(sat⁺, sens)

            # Take a step 
            x̂⁺, d̂⁺, P̂⁺ = Estimator.sqrt_mekf(sat.state, sat.diodes, sat.covariance, alb, sensors.gyro, truth.Bᴵ, truth.sᴵ, 
                                             sensors.magnetometer, sensors.diodes, sensors.pos, dt, 
                                             data.Wchol, data.Vchol; calibrate_diodes = true)

            nP⁻ = norm(sat.covariance[7:end, 7:end])
            nP⁺ = norm(P̂⁺[7:end, 7:end])

            @test nP⁺ < nP⁻
        end
    end
end

#######################
# TEST FRAME OF OMEGA #
#######################
"""
    @testset "Omega frame" begin
        # x₀, _ = get_initial_state(; bias_less = true)
        # x₀ = STATE(x₀.r, x₀.v, 
        #             SVector{4, Float64}([cos(deg2rad(45)); [0, 1, 0] * sin(deg2rad(45))] ), 
        #             # SVector{4, Float64}([1, 0, 0, 0]),
        #             SVector{3, Float64}(deg2rad(0.0), deg2rad(10.0), deg2rad(0.0)), x₀.β)
        x₀ = STATE(
            [3e6; 3.8e6; -4.8e6],   # r 
            [-6861.0; 2331; -2476],   # v 
            [cos(deg2rad(45)); [0, 1, 0] * sin(deg2rad(45))],   # q 
            [deg2rad(10.0); deg2rad(0.0); deg2rad(0.0)],        # ω
            zeros(3)
        )
        
        vᴵ = [1, 0.0, 0]
        J  = SMatrix{3, 3, Float64, 9}([1.0 0 0; 0 1 0; 0 0 1])
        t  = Epoch(2020, 1, 1)
        u  = SVector{3, Float64}(0.0, 0.0, 0.0)
        dt = 0.1

        sat_state₀ = SAT_STATE(x₀.q, x₀.β)

        N = 360 
        vᴮ = zeros(N, 3);
        v̂ᴮ = zeros(N, 3);
        qs = zeros(N, 4);
        q̂s = zeros(N, 4);
        ωs = zeros(N, 3);
        x = x₀ 
        ss = sat_state₀
        for i = 1:N
            vᴮ[i, :] .= (H' * L(x.q) * R(x.q)' * H)' * vᴵ;
            qs[i, :] .= x.q
            global ss, _ = Estimator.prediction(ss, x.ω, dt, 6; calibrate_diodes = false)
            v̂ᴮ[i, :] .= (H' * L(ss.q) * R(ss.q)' * H)' * vᴵ;
            q̂s[i, :] .= ss.q
            ωs[i, :] .= x.ω
            global x = rk4(J, x, u, t, dt; σβ = 0.0);
        end

        # plot( [0, 1], [0, 0], [0, 0], label = "+X", c = :red,   lw = 3);
        # plot!([0, 0], [0, 1], [0, 0], label = "+Y", c = :blue,  lw = 3);  
        # plot!([0, 0], [0, 0], [0, 1], label = "+Z", c = :green, lw = 3, 
        #             xlim = [-1.1, 1.1], ylim = [-1.1, 1.1], zlim = [-1.1, 1.1] );
        # @gif for i ∈ 1:N
        #     plot!([0, vᴮ[i, 1]], [0, vᴮ[i, 2]], [0, vᴮ[i, 3]], c = :violet, lw = 2, label = false)
        #                 # xlim = [-1.1, 1.1], ylim = [-1.1, 1.1], zlim = [-1.1, 1.1] )
        # end every 5

        # p1 = plot(vᴮ, title = "RK4", label = ["x" "y" "z"], c = [:red :blue :green])
        # p2 = plot(v̂ᴮ, title = "Pred", label = ["x" "y" "z"], c = [:red :blue :green])
        # plot(p1, p2, layout = 2)
    end
"""


