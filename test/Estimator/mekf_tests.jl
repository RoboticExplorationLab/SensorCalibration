# [test/Estimator/mekf_tests.jl]

""" To Do:
 - A (in prediction) does not appear to be correct
 - H in mag_meas is half? Same with current_meas
 - These are all kinda messy tests (pull out shared things, delete useless)

 - I am a little concerned that the covariance doesnt decrease in the last of 'estimate'

 - in mekf, we can't use x.ω, we need x.ω + β
 - Check for ===
"""


######## Stuff needed to run this independently ##################################
# using Test, BenchmarkTools 
# using StaticArrays, LinearAlgebra
# using ForwardDiff, Plots, EarthAlbedo, SatelliteDynamics, JLD2

# include("../../src/CustomStructs.jl");   using .CustomStructs 
# include("../../src/Estimator/mekf.jl");

# include("../../src/quaternions.jl")
# include("../../src/mag_field.jl")
# # # (Also change the path in get_albedo)
##################################################################################

@testset "MEKF Tests" begin 
    import .Estimator.prediction, .Estimator.mag_measurement, .Estimator.current_measurement
    import .Estimator.sqrt_mekf, .Estimator.estimate, .Estimator.reset_cov!

    function get_albedo(scale = 1) 

        function load_refl(path = "data/refl.jld2", scale = 1)
            temp = load(path)
        
            refl = REFL( temp["data"][1:scale:end, 1:scale:end], temp["type"], temp["start_time"], temp["stop_time"])
        
            return refl
        end
        lat_step = 1.0 * scale
        lon_step = 1.25 * scale

        refl = load_refl("../src/data/refl.jld2", scale)  
        cell_centers_ecef = get_albedo_cell_centers(lat_step, lon_step) 
        return ALBEDO(refl, cell_centers_ecef)
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

    function random_pose()
        p = randn(3)
        p = p / maximum(abs.(p))
        p =  6.871e6  * SVector{3, Float64}(p) #
    end;

    function vectrify(x::SAT_STATE{T}, d::DIODES{N, T}) where {N, T}
        v = zeros(T, 25)
        v[1:4]    .= x.q
        v[5:7]    .= x.β 
        v[8:13]   .= d.calib_values
        v[14:19]  .= d.azi_angles
        v[20:25]  .= d.elev_angles
        return v
    end;

    
    # Everything works except Jacobian A 
    @testset "Prediction" begin 

        # Standard function call (without diode)
        x = SAT_STATE()
        ω  = SVector{3, Float64}(randn(3))
        dt = 1.0
        N = 6
        x⁺, A =  prediction(x, ω, dt, N; calibrate_diodes = false)
        @test size(A) == (6, 6)
        @test size(x⁺.q) == (4,)
        @test norm(x⁺.q) ≈ 1.0
        # @code_warntype prediction(x, ω, dt)
        # @btime prediction($x, $ω, $dt; calibrate_diodes = true)


        # Standard function call, with diodes (shouldn't really be different, except size of A)
        x⁺_uncal = deepcopy(x⁺)   # To compare to 
        x2⁺, A =  prediction(x, ω, dt, N; calibrate_diodes = true)
        @test size(A) == (24, 24)
        @test all(A[7:end, 1:6] .== 0)
        @test all(A[1:6, 7:end] .== 0)
        @test A[7:end, 7:end] == I(18)
        @test (x⁺_uncal == x⁺) #&& (x⁺_uncal !== x⁺)     

        # No angular velocity or bias results in no change in q (without diodes)
        x = SAT_STATE(; β = zeros(3))
        ω  = SVector{3, Float64}(zeros(3))
        dt = 1.2
        x⁺, A =  prediction(x, ω, dt, N; calibrate_diodes = false)
        @test x.q == x⁺.q
        @test sum(diag(A)) ≈ 6
        
        # ... with diode calibration
        x = SAT_STATE(; β = zeros(3))
        ω  = SVector{3, Float64}(zeros(3))
        dt = 0.67
        x⁺, A =  prediction(x, ω, dt, N; calibrate_diodes = true)
        @test x.q == x⁺.q
        @test sum(diag(A)) == 24


        ##### Test state Jacobian A = ∂f/∂x #####
        #    Test that Q is SO(3), matches the quaternion provided; I, zeros in right spot
        x = SAT_STATE()
        ω  = SVector{3, Float64}(randn(3))
        dt = 0.92
        x⁺, A =  prediction(x, ω, dt, N; calibrate_diodes = false)
        Q = A[1:3, 1:3]
        @test (Q' * Q) ≈ I(3)
        @test det(Q) ≈ 1.0
        @test A[4:6, 4:6] == I(3)
        @test all(A[4:6, 1:3] .== 0.0)

        x = SAT_STATE()
        ω  = SVector{3, Float64}(randn(3))
        dt = 5.0
        x⁺, A =  prediction(x, ω, dt, N; calibrate_diodes = true)
        Q = A[1:3, 1:3]
        @test (Q' * Q) ≈ I(3)
        @test det(Q) ≈ 1.0
        @test A[4:end, 4:end] == I(21)
        @test all(A[4:6, 1:3] .== 0.0)


        #    Does it match ForwardDiff? <- No. See testset below 



        
        ##### Test x⁺ #####
        #    q is changed, but nothing else
        x = SAT_STATE()
        ω  = SVector{3, Float64}(randn(3))
        dt = 0.92
        x⁺, A =  prediction(x, ω, dt, N) 
        @test (x⁺.β == x.β) &&  x⁺.q ≉ x.q  
    end;

    # Our H is 0.5 * ForwardDiff's H
    @testset "Mag Measurement" begin 
        # Standard function call (without calibration)
        N = 6
        x  = SAT_STATE()
        Bᴵ = SVector{3, Float64}(randn(3))
        Bᴮ, Hb =  mag_measurement(x, Bᴵ, N; calibrate_diodes = false)
        # @code_warntype  mag_measurement(x, Bᴵ, N)
        # @btime  mag_measurement($x, $Bᴵ)         # Most of cost comes from quaternion functions
        @test (size(Bᴮ) == (3,)) && (size(Hb) == (3, 6))

        # Standard function call (with calibration) -> Should only affect sizes 
        x  = SAT_STATE()
        Bᴵ = SVector{3, Float64}(randn(3))
        Bᴮ, Hb =  mag_measurement(x, Bᴵ, N; calibrate_diodes = true)
        @test (size(Bᴮ) == (3,)) && (size(Hb) == (3, 24))


        ##### Test Bᴮ ##### 
        #  Bᴵ == Bᴮ if no rotation
        x  = SAT_STATE(; q = [1, 0, 0, 0])  # No rotation
        Bᴵ = SVector{3, Float64}(randn(3))
        Bᴮ, _ =  mag_measurement(x, Bᴵ, N)
        @test Bᴮ == Bᴵ 
        
        #  ||Bᴵ|| = ||Bᴮ|| even with rotation
        Nₜ = 100
        ts = zeros(Nₜ)
        for i = 1:Nₜ
            x  = SAT_STATE()  
            Bᴵ = SVector{3, Float64}(randn(3))
            Bᴮ, Hb =  mag_measurement(x, Bᴵ, N)
            ts[i] = (norm(Bᴮ) ≈ norm(Bᴵ)) && (quat2rot(x.q) * Bᴮ ≈ Bᴵ)
        end 
        @test sum(ts) == Nₜ

        #  Calibrating doesn't affect B 
        x  = SAT_STATE()  
        Bᴵ = SVector{3, Float64}(randn(3))
        Bᴮ_nocal, _ =  mag_measurement(x, Bᴵ, N; calibrate_diodes = false)
        Bᴮ_cal, _   =  mag_measurement(x, Bᴵ, N; calibrate_diodes = true)
        @test Bᴮ_nocal == Bᴮ_cal


        ##### Test H #####

        #   Test that it is [skew-symmetric  zeros] for no calibration
        Nₜ = 100
        ts = zeros(Nₜ)
        for i = 1:Nₜ
            x  = SAT_STATE()  # with rotation
            Bᴵ = SVector{3, Float64}(randn(3))
            Bᴮ, Hb =  mag_measurement(x, Bᴵ, N; calibrate_diodes = false)
            ts[i] = all(Hb[:, 4:6] .== 0.0) && all(diag(Hb[:, 1:3]) .== 0.0) &&
                        (Hb[1, 2] == -Hb[2, 1]) && (Hb[1, 3] == -Hb[3, 1])
        end 
        @test sum(ts) == Nₜ

        #   Test that it is [skew-symmetric  zeros] for calibration (just bigger zeros)
        Nₜ = 100
        ts = zeros(Nₜ)
        for i = 1:Nₜ
            x  = SAT_STATE()  # with rotation
            Bᴵ = SVector{3, Float64}(randn(3))
            Bᴮ, Hb =  mag_measurement(x, Bᴵ, N; calibrate_diodes = false)
            ts[i] = all(Hb[:, 4:end] .== 0.0) && all(diag(Hb[:, 1:3]) .== 0.0) &&
                        (Hb[1, 2] == -Hb[2, 1]) && (Hb[1, 3] == -Hb[3, 1])
        end 
        @test sum(ts) == Nₜ


        ### Test that it matches ForwardDiff ###
        # (1a) Create an FD-able version (remove struct, add in eltype, etc...)
        function mag_measurement_vec(x::Vector, Bᴵ::SVector{3, T}; N = 6, calibrate_diodes::Bool = true) where {T}
            q = x[1:4]
            ## Generate expected measurements
            Bᴮ = quat2rot(q)' * Bᴵ;  # Transpose to convert from inertial → body
        
            ## Generate Jacobian H (mag_measurement wrt state) 
        
            B̂ᴮ = hat(Bᴮ)  # Skew-symmetric matrix 
        
            H = (calibrate_diodes) ? zeros(eltype(x), 3, 6 + 3 * N) :  # Jacobian Matrix  H = [∂θ ∂β ∂C ∂α ∂ϵ]
                                     zeros(eltype(x), 3, 6)            # Jacobian Matrix  H = [∂θ ∂β] 
        
            H[1:3, 1:3] .= B̂ᴮ    
            # None of the other states directly affect measured mag field and are left as zeros
        
            return Bᴮ, H
        end;


        # (1b) Ensure this version matches ours
        Nₜ = 100 
        ts = zeros(Nₜ)
        for i = 1:Nₜ
            x = SAT_STATE()  # with rotation
            d = DIODES()
            Bᴵ = SVector{3, Float64}(randn(3))
            Bᴮ₁, H₁ =  mag_measurement(x, Bᴵ, N; calibrate_diodes = true)
            Bᴮ₂, H₂ = mag_measurement_vec(vectrify(x, d)[1:7], Bᴵ; calibrate_diodes = true)
            ts[i] = (Bᴮ₁ ≈ Bᴮ₂) && (H₁ ≈ H₂)
        end
        @test sum(ts) == Nₜ

        # (2) Use ForwardDiff and compare it to ours 
        x = SAT_STATE() 
        d = DIODES()
        Bᴵ = SVector{3, Float64}(randn(3))

        _mag_measurement_vec(_x) = mag_measurement_vec(_x, Bᴵ; calibrate_diodes = true)[1];  # Only want Bᴮ, not H
        Bᴮ_cust, H_cust =  mag_measurement(x, Bᴵ, N; calibrate_diodes = true)
        H_diff = ForwardDiff.jacobian(_mag_measurement_vec, vectrify(x, d))
        E(q) = [G(q)   zeros(4, 21);
                zeros(21, 3)   I(21) ];

        H_diff_adj = H_diff * E(x.q)
        
        @test H_diff_adj ≈ H_cust
        @test H_diff_adj ≈ 2 * H_cust
    end;

    # Our H is also 0.5 * ForwardDiff (for part)
    @testset "Current Measurement" begin 
        # Standard function call (without calibration)
        x = SAT_STATE() 
        d = DIODES()
        t = Epoch(2020, 4, 4)
        sᴵ = SVector{3, Float64}(sun_position(t))
        pos = 6.871e6 * SVector{3, Float64}([1.0; randn(2)])
        alb = get_albedo(1)
        Is, Hi =  current_measurement(x, d, sᴵ, pos, alb; calibrate_diodes = false)
        # @code_warntype  current_measurement(x, sᴵ, pos, alb)
        # @btime  current_measurement($x, $xᴵ, $pos, $alb);  # Slow, on order of ~most
        @test (size(Is) == (6,)) && (size(Hi) == (6, 6))
        @test all(Hi[Is .≤ 0.0, :] .== 0.0)

        # Standard function call (without albedo)
        x = SAT_STATE() 
        d = DIODES()
        t = Epoch(2020, 10, 10)
        sᴵ = SVector{3, Float64}(sun_position(t))
        pos = 6.871e6 * SVector{3, Float64}([randn(); 1.0; randn()])
        Is, Hi =  current_measurement(x, d, sᴵ, pos, alb; use_albedo = false)
        @test (size(Is) == (6,)) && (size(Hi) == (6, 24))
        @test all(Hi[Is .≤ 0.0, :] .== 0.0)

        # Standard function call (with albedo and calibration)
        x = SAT_STATE() 
        d = DIODES()
        t = Epoch(2020, 12, 25)
        sᴵ = SVector{3, Float64}(sun_position(t))
        pos = 6.871e6 * SVector{3, Float64}([randn(); randn(); 1.0])
        Is, Hi =  current_measurement(x, d, sᴵ, pos, alb; use_albedo = false)
        @test (size(Is) == (6,)) && (size(Hi) == (6, 24))
        @test all(Hi[Is .≤ 0.0, :] .== 0.0)
    
 
        ##### Test Is (diode currents) #####
        #    Non-negative, both with and without albedo 
        x = SAT_STATE() 
        d = DIODES()
        t = Epoch(2021, 7, 1)
        sᴵ = SVector{3, Float64}(sun_position(t))
        pos = random_pose()
        alb = get_albedo(1)
        Is, _ =  current_measurement(x, d, sᴵ, pos, alb; use_albedo = true)
        @test all(Is .≥ 0.0)

        x = SAT_STATE() 
        d = DIODES()
        t = Epoch(2019, 10, 15)
        sᴵ = SVector{3, Float64}(sun_position(t))
        pos = random_pose()  
        Is, _ =  current_measurement(x, d, sᴵ, pos, alb; use_albedo = false)
        @test all(Is .≥ 0.0)

        #    One per pair has to be zero for evenly spaced sides and NO albedo
        x = SAT_STATE(; ideal = true);
        d = DIODES(; ideal = true);
        t = Epoch(2019, 7, 1);
        sᴵ = SVector{3, Float64}(sun_position(t));
        pos = random_pose();
        Is, _ =  current_measurement(x, d, sᴵ, pos, alb; use_albedo = false);
        @test length(Is[Is .> 0.001]) == 3 
        
        ##### Test measurement Jacobian H #####
        #    H is 0 where I is 0; ∂β is zero
        x = SAT_STATE() 
        d = DIODES()
        t = Epoch(2020, 5, 4)
        sᴵ = SVector{3, Float64}(sun_position(t))
        pos = random_pose()
        alb = get_albedo(2)
        Is, Hi =  current_measurement(x, d, sᴵ, pos, alb)
        @test all(Hi[:, 4:6] .== 0.0)
        @test all(Hi[Is .== 0.0, :] .== 0.0)

        x = SAT_STATE() 
        d = DIODES()
        t = Epoch(2020, 5, 4)
        sᴵ = SVector{3, Float64}(sun_position(t))
        random_pose
        Is, Hi =  current_measurement(x, d, sᴵ, pos, alb; use_albedo = false)
        @test all(Hi[:, 4:6] .== 0.0)
        @test all(Hi[Is .== 0.0, :] .== 0.0)


        
        ##  Comparison to old version (Need a few helper functions)
        function current_measurement_old(x, 𝐬ᴵ, i, pos, time, alb::ALBEDO; _E_am0 = 1366.9)  
        
            q, β, c, α, ϵ = split_state(x, i)

            ᴮQᴵ = dcm_from_q(q)'; # DCM from quaternion (transposed to get I -> B)
            sᴮ = ᴮQᴵ * 𝐬ᴵ 
        
            šᴮ= hat(sᴮ);  # Skew-symmetric form
            n   = [cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)];  # [i x 3]
            ndα = [(-cos.(ϵ).*sin.(α)) (cos.(ϵ).*cos.(α)) zeros(size(α))];
            ndϵ = [(-sin.(ϵ).*cos.(α)) ((-sin.(ϵ).*sin.(α))) cos.(ϵ)]; # (With negative middle term, differing from the paper)
            
            ∂θ = (c .* n) * šᴮ;     # [i x 3]
            ∂β = zeros(i, 3);       # [i x 3]  
            ∂C = n * sᴮ;            # [i,]
            ∂α = c .* (ndα * sᴮ);   # [i,]
            ∂ϵ = c .* (ndϵ * sᴮ);   # [i,]  
        
            H = [∂θ ∂β Diagonal(∂C) Diagonal(∂α) Diagonal(∂ϵ)] # [i x 6 + 3i]
            I_meas = c .* (n * sᴮ)  # Measured current, ALBEDO added in later
        
        
            # ADD IN ALBEDO
            sᴵ_unscaled = sun_position(time) - pos;
            # ecl = eclipse_conical(-pos, sᴵ_unscaled) ####### NEED TO FIX TO +pos when updated
            # ecl = (ecl > 0.98) ? 1.0 : 0.0
        
            sPos = SVector{3, Float64}(pos)
            ssᴵ  = SVector{3, Float64}(sᴵ_unscaled)
            albedo_matrix = earth_albedo(sPos, ssᴵ, alb.refl.data) 
        
        
        
            for j = 1:i
                # surface_normal = [cos(ϵ[j])*cos(α[j]) cos(ϵ[j])*sin(α[j]) sin(ϵ[j])]     # Photodiode surface normal
                
                diode_albedo = compute_diode_albedo_old(albedo_matrix, alb.cell_centers_ecef, view(n, j, :), pos)
                 
                diode_albedo = c[j] * diode_albedo / _E_am0;
                I_meas[j] = I_meas[j] + diode_albedo
            end
            
            # Account for eclipses
            # I_meas *= ecl
            I_meas[I_meas .< 0] .= 0  # Photodiodes don't generate negative current
            H[I_meas .≤ 0, :] .= 0    # ^ To match the above
            y = I_meas[:]             # [i x 1]        
            return y, H
        end

        function split_state(x, i)
            """ Helper function to split state into various components """
            x = x[:]
            q = x[1:4]
            β = x[5:7]
            c = x[8:(7+i)]  
            α = x[(8+i):(7+2*i)]
            ϵ = x[(8+2*i):end]
        
            return q, β, c, α, ϵ
        end

        function compute_diode_albedo_old(albedo_matrix, cell_centers_ecef, surface_normal, sat_pos)
            """ 
                Estimates the effect of Earth's albedo on a specific photodiode (by using the surface normal of that diode)
                    = cell_albedo * surface_normal^T * r_g, with r_g as a unit vector in the direction of the grid point on Earth
        
                Arguments:
                - albedo_matrix: Albedo values for each cell on the Earth's surface         | [num_lat x num_lon] 
                - surface_normal: Photodiode surface normal                                 | [3,]
                - sat_pos: Cartesian position of satellite                                  | [3,]
        
                Returns:
                - diode_albedo: Total effect of albedo on specified photodiode              | Scalar
            """    
            diode_albedo = 0.0
            r_g = zeros(Float64, 3)
            for r = 1:1:size(albedo_matrix, 1)
                for c = 1:1:size(albedo_matrix, 2)
                    if albedo_matrix[r,c] != 0
                        r_g .= view(cell_centers_ecef, r, c, :) .- sat_pos
                        r_g .= r_g ./ norm(r_g)  # Make unit
        
                        cell_albedo = (albedo_matrix[r,c] * dot(surface_normal, r_g))
        
                        if cell_albedo > 0.0    # Can't be negative
                            diode_albedo += cell_albedo 
                        end
                    end
                end
            end
            
            return diode_albedo
        end

        function dcm_from_q(quat)
            # Takes in a quaternion (scalar last - [q⃗, q]) and returns the DCM  (From Kevin)
            
            quat = quat[:];
            
            v = quat[1:3];
            s = quat[4];
        
            DCM = I(3) + 2*hat(v)*(s*I(3) + hat(v));
        
            return DCM;
        end


        # Compare old to new 
        Nt = 100
        ts = zeros(Nt)
        for i = 1:Nt
            x = SAT_STATE() 
            d = DIODES()
            t = Epoch(2020, 5, 4)
            sᴵ = SVector{3, Float64}(sun_position(t))
            pos = random_pose()
            alb = get_albedo(1)

            Is_new, H_new =  current_measurement(x, d, sᴵ, pos, alb; calibrate_diodes = true, use_albedo = true)

            xVec = vectrify(x, d)
            xVec[1:4] .= [xVec[2:4]; xVec[1]]
            Is_old, H_old = current_measurement_old(xVec, (sᴵ / norm(sᴵ)), 6, pos, t, alb)

            ts[i] = (Is_new ≈ (Is_old / norm(Is_old))) && (H_new  ≈ H_old)
        end
        @test all(ts .== 1)


        ### Test that it matches ForwardDiff ### 
        # (1a) Create an FD-able version (remove struct, add eltype(x))
        function current_measurement_vec(x::Vector, sᴵ::SVector{3, T}, pos::SVector{3, T}, alb::ALBEDO; 
            E_am₀ = 1366.9, use_albedo = true, calibrate_diodes::Bool = true, N = 6) where {T}
            
            q, β, C, α, ϵ = split_state(x, N)

            # Generate expected measurements
            sᴵ_unit = sᴵ / norm(sᴵ)  # Unit sun vector in inertial frame
            sᴮ = quat2rot(q)' * (sᴵ_unit);  # Transpose to convert from inertial → body
            n =([cos.(ϵ).*cos.(α) cos.(ϵ).*sin.(α) sin.(ϵ)]);  # Surface normals for each photodiode
            currents = C .* (n * sᴮ)     # Expected current measurements

            albedos  = zeros(eltype(x), N)
            if use_albedo
                albedo_matrix = earth_albedo(pos, sᴵ, alb.refl.data)
                for i = 1:N 
                    diode_albedo = compute_diode_albedo(albedo_matrix, alb.cell_centers_ecef,  n[i, :] , pos)  
                    albedos[i] = (C[i] / E_am₀) * diode_albedo 
                end
            end

            # Account for eclipses in I, H (a bit problematic for H)
            Is = zeros(eltype(x), N)  # This is to deal with the pesky Static problems
            for i = 1:N 
                Is[i] = currents[i] + albedos[i]
            end

            # # Generate Jacobian H (current measurement wrt state)
            ŝᴮ = hat(sᴮ);          # Skew-symmetric 
            ∂θ = (C .* n) * ŝᴮ   # How slight changes in orientation affect expected current field (∂I/∂θ)
            ∂β = zeros(T, N, 3)    # How slight changes in gyro bias affect expected current field (∂I/∂β)


            # H = [∂θ ∂β ∂C ∂α ∂ϵ]  :  [∂θ ∂β]
            H = (calibrate_diodes) ? zeros(eltype(x), 6, 6 + 3 * N)  : zeros(eltype(x)s, 6, 6)  
            H[:, 1:3] .= ∂θ
            H[:, 4:6] .= ∂β    # Should be just zeros...


            if calibrate_diodes
                # Jacobian Matrix  H = [∂θ ∂β ∂C ∂α ∂ϵ]
                ndα = [(-cos.(ϵ).*sin.(α)) ( cos.(ϵ).*cos.(α)) zeros(size(α))];
                ndϵ = [(-sin.(ϵ).*cos.(α)) (-sin.(ϵ).*sin.(α)) cos.(ϵ)]; # (With negative middle term, differing from the paper)

                ∂C = n * sᴮ;            # [i,]
                ∂α = C .* (ndα * sᴮ);   # [i,]
                ∂ϵ = C .* (ndϵ * sᴮ);   # [i,]  

                H[:, 7:end] .= [Diagonal(∂C) Diagonal(∂α) Diagonal(∂ϵ)]
            end

            Is[Is .< 1e-8]   .= 0.0   # Photodiodes don't generate negative current
            H[Is .≤  0.0, :] .= 0.0   #   ^ to match the above 

            return Is, H
        end;

        function compute_diode_albedo(data::Matrix{T}, cell_centers_ecef::Array{T, 3}, surface_normal, sat_pos::SVector{3, T}) where {T}

            Nlat, Nlon = size(data)
        
            diode_albedo = 0.0
            r_g = zeros(eltype(surface_normal), 3)
        
            for r = 1:Nlat
                for c = 1:Nlon
        
                    if data[r,c] != 0
                        r_g .= view(cell_centers_ecef, r, c, :) .- sat_pos
                        r_g .= r_g ./ norm(r_g)  # Make unit
        
                        cell_albedo = (data[r,c] * dot(surface_normal, r_g))
        
                        if cell_albedo > 0.0    # Can't be negative
                            diode_albedo += cell_albedo 
                        end
                    end
                end
            end
            
            return diode_albedo
        end;


        # (1b) Compare this to ours
        Nₜ = 50
        ts = zeros(Nₜ)
        alb = get_albedo(2)
        for i = 1:Nₜ
            x = SAT_STATE() 
            d = DIODES()
            t = Epoch(2020, 5, 4)
            sᴵ = SVector{3, Float64}(sun_position(t))
            pos = random_pose()
            Is, Hi =  current_measurement(x, d, sᴵ, pos, alb; calibrate_diodes = true)
            xVec = vectrify(x, d)
            Is2, Hi2 = current_measurement_vec(xVec, sᴵ, pos, alb)
            ts[i] = (Is ≈ Is2) && (Hi ≈ Hi2)
        end 
        @test all(ts .== 1.0)


        # (2) Compare to ForwardDiff (note custom jacobian ignores albedo)
        _current_measurement_vec(_x) = current_measurement_vec(_x, sᴵ, pos, alb; calibrate_diodes = true, use_albedo = false)[1];  # Only want Is

        x = SAT_STATE() 
        d = DIODES()
        t = Epoch(2020, 5, 4)
        sᴵ = SVector{3, Float64}(sun_position(t))
        pos = random_pose()
        Is, Hi =  current_measurement(x, d, sᴵ, pos, alb; calibrate_diodes = true, use_albedo = false)
        H_diff = ForwardDiff.jacobian(_current_measurement_vec, vectrify(x, d))
       
        E(q) = [G(q)   zeros(4, 21);
                zeros(21, 3)   I(21) ];

        # The calibration part matches but the rest is 1/2 off
        H_diff_adj = H_diff * E(x.q)
        @test H_diff_adj[:, 1:3] ≈ Hi[:, 1:3]
        @test H_diff_adj[:, 1:3] ≈ 2 * Hi[:, 1:3]

        @test H_diff_adj[:, 4:end] ≈ Hi[:, 4:end]
    end;

    # Not really sure how to test this
    @testset "Sqrt MEKF" begin 

        # LOTS of data to prep
        x = SAT_STATE(; ideal = true)
        d = DIODES(; ideal = true)   # NOTE this 'ideal' assumes the 45* rotation so ±Z isnt unobservable
        σq = deg2rad(10); σβ = deg2rad(10)
        Σϕ = diagm( σq * ones(3) )  
        Σβ = diagm( σβ * ones(3) )
        Σ  = I(6)
        Pchol = SAT_COVARIANCE(Σϕ, Σβ, Σ, Σ, Σ).Σ[1:6, 1:6]

        alb = get_albedo(2)
        Bᴵ  = SVector{3, Float64}(randn(3))
        sᴵ  = SVector{3, Float64}(AU * [1.0, 0.0, 0.0])
        B̃ᴮ  = Bᴵ 
        Ĩ   = SVector{6, Float64}(sqrt(2)/2, 0, 0, 0, sqrt(2)/2, 0) # Should be all zeros because of eclipse, but that isnt checked for 
        pos = SVector{3, Float64}(6.871e6 * [-1.0, 0.0, 0.0])
        dt = 1.0
        ω   = SVector{3, Float64}(zeros(3))
        N = 6

        gyro_bias_instability = 0.8
        angle_random_walk     = 0.06
        σ_gyro = deg2rad(gyro_bias_instability) / 3600.0  # Convert (deg/hour) to (rad/sec)
        σ_bias = deg2rad(angle_random_walk) / 60.0        # Convert (deg/sqrt(hour)) to ( rad/sqrt(s) )
        σ_sun = deg2rad(3.0);
        σ_mag = deg2rad(3.0);

        W = Diagonal( [σ_gyro * ones(3); σ_bias * ones(3)].^2 )
        Wchol = chol(Matrix(W)) 
        V = Diagonal( [σ_mag * ones(3); σ_sun * ones(N)].^2 )
        Vchol = chol(Matrix(V))

        x⁺, d⁺, Pchol⁺ =  sqrt_mekf(x, d, Pchol, alb, ω, Bᴵ, sᴵ, B̃ᴮ, Ĩ, pos, dt, Wchol, Vchol; calibrate_diodes = false);

        # No spin, so quat should be the same 
        @test x⁺.q ≈ x.q  atol = 1e-14
        @test x⁺.β ≈ x.β  atol = 1e-14
        @test d⁺.calib_values == d.calib_values   # These three should be untouched
        @test d⁺.azi_angles == d.azi_angles 
        @test d⁺.elev_angles == d.elev_angles

        @test norm(Pchol) > norm(Pchol⁺)  # Because our measurements and dynamics match
    end

    @testset "Estimate" begin 
        x = SAT_STATE(; ideal = true )
        d = DIODES(; ideal = true)  # Assumes `ideal` is with 45* rotation (affects Ĩ )
        sat = SATELLITE(; sta = x, dio = d)
        alb = get_albedo(2) 
        t   = Epoch(2021, 5, 2)
        dt  = 1.0

        sᴵ  = SVector{3, Float64}(sun_position(t))
        pos = SVector{3, Float64}(6.871e6 * [-1.0, 0.0, 0.0])
        B̃ᴮ  = SVector{3, Float64}(IGRF13(pos, t))
        Ĩ   = SVector{6, Float64}(0.3379, 0.0, 0.61329, 0.0, 0.71393, 0.0) # Computed beforehand 
        dt = 1.0
        ω   = SVector{3, Float64}(zeros(3))
        N = 6

        gyro_bias_instability = 0.8
        angle_random_walk     = 0.06
        σ_gyro = deg2rad(gyro_bias_instability) / 3600.0  # Convert (deg/hour) to (rad/sec)
        σ_bias = deg2rad(angle_random_walk) / 60.0        # Convert (deg/sqrt(hour)) to ( rad/sqrt(s) )
        σ_sun = deg2rad(3.0);
        σ_mag = deg2rad(3.0);

        W = Diagonal( [σ_gyro * ones(3); σ_bias * ones(3)].^2 )
        Wchol = chol(Matrix(W)) 
        V = Diagonal( [σ_mag * ones(3); σ_sun * ones(N)].^2 )
        Vchol = chol(Matrix(V))

        sens = SENSORS(B̃ᴮ, Ĩ, ω, pos)
        data = MEKF_DATA(Wchol, Vchol)

        sat⁺ =  estimate(sat, sens, data, alb, t, dt; calibrate_diodes = false);

        @test sat.covariance !== sat⁺.covariance   # Make sure we have the right copy
        @test sat.state   !==  sat⁺.state

        @test sat.J == sat⁺.J 
        @test sat.diodes == sat⁺.diodes 
        @test sat.magnetometer == sat⁺.magnetometer 

        cov = sat.covariance; cov⁺ = sat⁺.covariance 
        state = sat.state; state⁺ = sat⁺.state
        diodes = sat.diodes; diodes⁺ = sat⁺.diodes
        @test cov[1:6, 1:6]  ≉ cov⁺[1:6, 1:6]
        @test cov[7:end, 7:end] ≈ cov⁺[7:end, 7:end]
        @test norm(cov) > norm(cov⁺)

        @test diodes.calib_values == diodes⁺.calib_values 
        @test diodes.azi_angles == diodes⁺.azi_angles
        @test diodes.elev_angles == diodes⁺.elev_angles 
        @test state.q ≈ state⁺.q  atol = 1e-5 # Because guess was right
        @test state.β ≈ state⁺.β  atol = 1e-5 # Because guess was right 



        #### TRIAL 2, with slightly not ideal
        q = [1.1, 0.05, 0.1, 0.05]; q = SVector{4, Float64}(q / norm(q))
        β = SVector{3, Float64}(0.1, -0.1, 0.05)
        x = SAT_STATE(; q = q, β = β)
        d = DIODES(; ideal = true) 
        sat = SATELLITE(; sta = x, dio = d)
        alb = get_albedo(2) 
        t   = Epoch(2021, 12, 2)
        dt  = 1.0

        sᴵ  = SVector{3, Float64}(sun_position(t))
        pos = SVector{3, Float64}(6.871e6 * [0.0, 0.0, -1.0])
        B̃ᴮ  = SVector{3, Float64}((IGRF13(pos, t)))
        Ĩ   = SVector{6, Float64}(0.1228, 0.2411, 0.0, 0.9812, 0.0, 0.5117) # Computed beforehand, with some noise
        dt = 2.0
        ω   = SVector{3, Float64}([0.1, -0.1, 0.2])
        N = 6

        sens = SENSORS(B̃ᴮ, Ĩ, ω, pos)
        data = MEKF_DATA()

        sat⁺ =  estimate(sat, sens, data, alb, t, dt; calibrate_diodes = false);

        @test sat.covariance !== sat⁺.covariance   # Make sure we have the right copy
        @test sat.state   !==  sat⁺.state

        @test sat.J == sat⁺.J 
        @test (sat.diodes == sat⁺.diodes) # && (sat.diodes !== sat⁺.diodes)
        @test sat.magnetometer == sat⁺.magnetometer 

        cov = sat.covariance; cov⁺ = sat⁺.covariance 
        state = sat.state; state⁺ = state⁺ = sat⁺.state
        diodes = sat.diodes; diodes⁺ = sat⁺.diodes
        @test cov[1:6, 1:6]  ≉ cov⁺[1:6, 1:6]
        @test cov[7:end, 7:end] ≈ cov⁺[7:end, 7:end]
        # @test norm(cov) > norm(cov⁺)

        @test diodes.calib_values == diodes⁺.calib_values 
        @test diodes.azi_angles == diodes⁺.azi_angles 
        @test diodes.elev_angles == diodes⁺.elev_angles 

        @test state.q ≉ state⁺.q  atol = 1e-6 # Because guess was not quite right
        @test state.β ≉ state⁺.β  atol = 1e-6 # Because guess was not quite right
        

        # TRIAL 3 with garbage values
        sat = SATELLITE()
        alb = get_albedo(2) 
        t   = Epoch(2021, 12, 2)
        dt  = 2.0

        sᴵ  = SVector{3, Float64}(sun_position(t))
        pos = SVector{3, Float64}(6.871e6 * [0.0, 0.0, -1.0])
        B̃ᴮ  = SVector{3, Float64}(randn(3))
        Ĩ   = SVector{6, Float64}(rand(6)) # Computed beforehand, with some noise
        ω   = SVector{3, Float64}([0.1, -0.1, 0.2])
        N = 6

        sens = SENSORS(B̃ᴮ, Ĩ, ω, pos)
        data = MEKF_DATA()

        sat⁺ =  estimate(sat, sens, data, alb, t, dt; calibrate_diodes = false);

        @test sat.covariance !== sat⁺.covariance   # Make sure we have the right copy
        @test sat.state   !==  sat⁺.state

        @test sat.J == sat⁺.J 
        @test sat.diodes == sat⁺.diodes 
        @test sat.magnetometer == sat⁺.magnetometer 

        cov = sat.covariance; cov⁺ = sat⁺.covariance 
        state = sat.state; state⁺ = sat⁺.state
        diodes = sat.diodes; diodes⁺ = sat⁺.diodes
        @test cov[7:end, 7:end] ≈ cov⁺[7:end, 7:end]
       
        @test diodes.calib_values == diodes⁺.calib_values 
        @test diodes.azi_angles == diodes⁺.azi_angles
        @test diodes.elev_angles == diodes⁺.elev_angles 




        # TRIAL 4 with garbage values and allowing for calibration
        sat = SATELLITE()
        alb = get_albedo(2) 
        t   = Epoch(2021, 12, 2)
        dt  = 2.0

        sᴵ  = SVector{3, Float64}(sun_position(t))
        pos = SVector{3, Float64}(6.871e6 * [0.0, 0.0, -1.0])
        B̃ᴮ  = SVector{3, Float64}(randn(3))
        Ĩ   = SVector{6, Float64}(rand(6)) # Computed beforehand, with some noise
        ω   = SVector{3, Float64}([0.1, -0.1, 0.2])
        N = 6
        
        gyro_bias_instability = 0.8
        angle_random_walk     = 0.06
        σ_gyro = deg2rad(gyro_bias_instability) / 3600.0  # Convert (deg/hour) to (rad/sec)
        σ_bias = deg2rad(angle_random_walk) / 60.0        # Convert (deg/sqrt(hour)) to ( rad/sqrt(s) )
        σ_mag = deg2rad(2.0);
        σ_cur = 0.1;

        W = Diagonal( [σ_gyro * ones(3); σ_bias * ones(3)].^2 )
        Wchol = chol(Matrix(W)) 
        V = Diagonal( [σ_mag * ones(3); σ_cur * ones(N)].^2 )
        Vchol = chol(Matrix(V))

        sens = SENSORS(B̃ᴮ, Ĩ, ω, pos)
        data = MEKF_DATA()
        @test data.Wchol[1:6, 1:6] == Wchol
        @test data.Vchol == Vchol
      
        sat⁺ =  estimate(sat, sens, data, alb, t, dt; calibrate_diodes = true);

        @test sat.covariance !== sat⁺.covariance   # Make sure we have the right copy
        @test sat.state   !==  sat⁺.state

        @test sat.J == sat⁺.J 
        @test sat.diodes != sat⁺.diodes 
        @test sat.magnetometer == sat⁺.magnetometer 

        cov = sat.covariance; cov⁺ = sat⁺.covariance 
        state = sat.state; state⁺ = sat⁺.state
        diodes = sat.diodes; diodes⁺ = sat⁺.diodes
        @test cov[7:end, 7:end] ≉ cov⁺[7:end, 7:end]
        
        @test diodes.calib_values != diodes⁺.calib_values 
        @test diodes.azi_angles != diodes⁺.azi_angles
        @test diodes.elev_angles != diodes⁺.elev_angles 





        # Test 5 - correct guess with diode cal
        x = SAT_STATE(; ideal = true )
        sat = SATELLITE(; sta = x)
        alb = get_albedo(2) 
        t   = Epoch(2021, 5, 2)
        dt  = 1.0

        sᴵ  = SVector{3, Float64}(sun_position(t))
        pos = SVector{3, Float64}(6.871e6 * [-1.0, 0.0, 0.0])
        B̃ᴮ  = SVector{3, Float64}(IGRF13(pos, t))
        Ĩ, _ = current_measurement(x, sat.diodes, sᴵ, pos, alb; use_albedo = true, calibrate_diodes = false)
        Ĩ   = SVector{6, Float64}(Ĩ)
        dt = 1.0
        ω   = SVector{3, Float64}(zeros(3))
        N = 6

        sens = SENSORS(B̃ᴮ, Ĩ, ω, pos)
        data = MEKF_DATA()

        sat⁺ =  estimate(sat, sens, data, alb, t, dt; calibrate_diodes = true)

        @test sat.covariance !== sat⁺.covariance   # Make sure we have the right copy
        @test sat.state   !==  sat⁺.state

        @test sat.J == sat⁺.J 
        # @test (sat.diodes == sat⁺.diodes) && (sat.diodes !== sat⁺.diodes) 
        @test sat.magnetometer == sat⁺.magnetometer 

        cov = sat.covariance; cov⁺ = sat⁺.covariance 
        state = sat.state; state⁺ = sat⁺.state
        diodes = sat.diodes; diodes⁺ = sat⁺.diodes

        @test cov[1:6, 1:6]  ≉ cov⁺[1:6, 1:6]
        @test cov[7:end, 7:end] ≉ cov⁺[7:end, 7:end]
        # @test norm(cov) > norm(cov⁺)

        @test diodes.calib_values ≈ diodes⁺.calib_values atol = 1e-3
        @test diodes.azi_angles ≈ diodes⁺.azi_angles atol = 1e-3 
        @test diodes.elev_angles ≈ diodes⁺.elev_angles atol = 1e-3 
        @test state.q ≈ state⁺.q atol = 1e-3   
        @test state.β ≈ state⁺.β atol = 1e-3    
    end

    @testset "Reset cov" begin 
        sat = SATELLITE()
        sat.covariance .= sat.covariance * sat.covariance
        cov₀ = deepcopy(sat.covariance)

         reset_cov!(sat; reset_calibration = false);
        @test cov₀ != sat.covariance

        sat = SATELLITE()
        sat.covariance .= sat.covariance * sat.covariance
        cov₀ = deepcopy(sat.covariance)

         reset_cov!(sat; reset_calibration = true);
        @test cov₀ != sat.covariance
    end
end



"""
    @testset "Attempt to figure out A in Pred..." begin 
        # FD-able version
        function prediction_FD(x::SAT_STATE{T}, ω::SVector{3, T}, dt::T, N::Int; calibrate_diodes::Bool = true) where {T}

            function nextq(q, ω, β)
                # Generate axis-angle representation
                γ  = ω - β     # Corrected angular velocity 
                nγ = norm(γ)     # Magnitude of corrected angular velocity 
            
                # Predict next orientation (as a quaternion)
                r  = γ / nγ       # Unit axis of rotation
                θ  = nγ * dt      # Angular rotaation about unit axis
                q⁺ = (nγ == 0.0) ?  q :
                                    qmult(q, [cos(θ / 2); r*sin(θ / 2)])

                return q⁺, γ, nγ
            end

            q⁺, γ, nγ = nextq(x.q, ω, x.β)

            x⁺ = update_state(x; q = q⁺)  # Doesn't matter if we are calibrating diodes, only q changes in prediction
        
            # Calculate Jacobian ∂f/∂x
        
            γ̂ = -hat(γ / nγ)  # Skew-symmetric matrix 
            R = (nγ == 0.0) ? I(3)  :  # Avoid the divide-by-zero error
                            I(3) + (γ̂ ) * sin(nγ) + ((γ̂ )^2) * (1 - cos(nγ)); # Rodrigues formula

            @debug R == exp(nγ * hat(γ / nγ) )
                            
            A = zeros(T, 6, 6)     # Jacobian of f(x) wrt x,   A  =   [R   -dt I₃;  0₃  I₃]


            _nextq_q(_q) = nextq(_q, ω, x.β)[1]
            _nextq_β(_β) = nextq(x.q, ω, _β)[1]

            q⁺ₐ = Array(q⁺) # To prevent FD from somehow overwriting q⁺...
            A[1:3, 1:3] .= attitude_jacobian(q⁺ₐ)' * ForwardDiff.jacobian(_nextq_q, x.q) * attitude_jacobian(x.q)   # df/dq with Attitude Jacobian 
            A[1:3, 4:6] .= attitude_jacobian(q⁺ₐ)' * ForwardDiff.jacobian(_nextq_β, x.β) 
            # A[1:3, 1:3] .= R 
            # A[1:3, 4:6] .= -dt * I(3) 
            A[4:6, 4:6] .= I(3)
            # A[4:6, 1:3] .= zeros(3, 3)
        
            
            if calibrate_diodes
                # Add in additional states if calibrating (i.e., x = [q, β, C, α, ϵ])
                Ac = zeros(T, 6 + 3 * N, 6 + 3 * N)  # Ac = [A  0; 0 I]
                Ac[1:6, 1:6] .= A 
                Ac[7:end, 7:end] .= I(3 * N)
                
                return x⁺, Ac
            else
                return x⁺, A
            end
        end;

        function attitude_jacobian(q)
            G = [-q[2:4]' ; 
                q[1] * I + hat(q[2:4])]
            return SMatrix{4, 3, Float64, 12}(G)
        end

        # MAKE SURE THE x⁺ MATCH
        x = SAT_STATE()
        ω = SVector{3, Float64}(randn(3))
        N = 6 
        dt = 0.1 

        xpred_fd, A_fd = prediction_FD(x, ω, dt, N; calibrate_diodes = false)
        xpred_og, A_og = prediction(x, ω, dt, N; calibrate_diodes = false)
        @test xpred_fd == xpred_og
        @test A_fd ≈ A_og

    end
"""
