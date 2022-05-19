# [test/Simulator/measurement_tests.jl]

""" To Do:
 - Verify that rotation_noise adds correct Gaussian
 - to avoid scope conflicts, I am just calling CustomStructs through Simulator rn...
 """

@testset "Measurement Tests" begin 

    function get_albedo(scale = 1) 

        function load_refl(path = "data/refl.jld2", scale = 1)
            temp = load(path)
        
            refl = REFL( temp["data"][1:scale:end, 1:scale:end], temp["type"], temp["start_time"], temp["stop_time"])
        
            return refl
        end
        lat_step = 1.0 * scale
        lon_step = 1.25 * scale

        refl = load_refl("../src/MissionSim/data/refl.jld2", scale)  # Not sure this will quite work...
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

    @testset "  Mag Calibration Matrix" begin 
        sat = Simulator.SATELLITE() 

        M = Simulator.get_mag_calibration_matrix(sat)
        @test size(M) == (3, 3)  

        # Lower triangular
        @test M[1, 2] == 0   
        @test M[1, 3] == 0 
        @test M[2, 3] == 0 

        @test M[1, 1] > 0.0
        @test M[1, 1] == sat.magnetometer.scale_factors[1]

        allocs = @ballocated Simulator.get_mag_calibration_matrix($sat) evals = 1 samples = 1
        @test allocs == 160 # bytes, size of a 3 x 3 matrix
    end;

    @testset "  Compute Diode Albedo" begin 
        function plot_albedo(cell_centers, vals, sat, sun)
            sun = sun / 1000
            sat = sat / 0.9

            l, w, N = size(cell_centers)
            cc = reshape(cell_centers, l * w, N)
            v  = reshape(vals, l * w)

            # Downsample 
            cc = cc[1:7:end, :]
            v  = v[1:7:end]

            lit = cc[v .> 0.01, :]
            unlit = cc[v .≤ 0.01, :]

            scatter(lit[:, 1], lit[:, 2], lit[:, 3], c = :yellow, markerstrokecolor = :yellow, ms = 1)
            scatter!(unlit[:, 1], unlit[:, 2], unlit[:, 3], ms = 1, c = :black)
            scatter!([sun[1]], [sun[2]], [sun[3]], c = :yellow, ms = 5)
            display(scatter!([sat[1]], [sat[2]], [sat[3]], c = :blue, ms = 2))

            return cc, v
        end

        alb = get_albedo()

        # Set up photodiode surface normals (ideally, we would have one on each side of the CubeSat)
        ns = [SVector{3, Float64}(1, 0, 0), SVector{3, Float64}(-1, 0, 0),
              SVector{3, Float64}(0, 1, 0), SVector{3, Float64}(0, -1, 0),
              SVector{3, Float64}(0, 0, 1), SVector{3, Float64}(0, 0, -1)]

        rsat = SVector{3, Float64}(6.5e6 * [0.5, 1.0, 0.1])
        rsun = SVector{3, Float64}(1e8 *   [1.0, 0.0, 0.1])

        data = earth_albedo(rsat, rsun, alb.refl.data) 
        cce  = alb.cell_centers_ecef
        # @code_warntype compute_diode_albedo(data, cce, ns[1], rsat)  
        alloc = @ballocated Simulator.compute_diode_albedo($data, $cce, $ns[1], $rsat);   # 100μs \ 1 alloc \ 112
        @test alloc == 112

        # # Simulator.Satellite is on opposite side
        
        # Trial 1
        rsun = SVector{3, Float64}(1e10  * [-1.0, 0.0, 0.0])
        rsat = SVector{3, Float64}(6.5e6 * [1.0, 0.0, 0.0])

        data = earth_albedo(rsat, rsun, alb.refl.data)   
        @test sum(data) == 0.0
        [@test Simulator.compute_diode_albedo(data, cce, ns[i], rsat) == 0.0 for i = 1:6]

        # Trial 2
        rsun = SVector{3, Float64}(9e9  * [0.0,  1.0, 0.0])
        rsat = SVector{3, Float64}(6.7e6 * [0.0, -1.0, 0.0])

        data = earth_albedo(rsat, rsun, alb.refl.data)   
        @test sum(data) == 0.0
        [@test Simulator.compute_diode_albedo(data, cce, ns[i], rsat) == 0.0 for i = 1:6]

        # Trial 3
        rsun = SVector{3, Float64}(1e11 * [ 0.5,  1.0,  0.3])
        rsat = SVector{3, Float64}(7e6  * [-0.5, -1.0, -0.3])

        data = earth_albedo(rsat, rsun, alb.refl.data)   
        @test sum(data) == 0.0
        [@test Simulator.compute_diode_albedo(data, cce, ns[i], rsat) == 0.0 for i = 1:6]

        # # Satellite is on same side

        # Trial 1 
        rsun = SVector{3, Float64}(1e10  * [1.0, 0.0, 0.0])
        rsat = SVector{3, Float64}(6.5e6 * [1.0, 0.0, 0.0])

        data = earth_albedo(rsat, rsun, alb.refl.data)  
        albs = [Simulator.compute_diode_albedo(data, cce, ns[i], rsat) for i = 1:6]
        @test albs[1] == 0.0            # Diode is facing the sun, not Earth
        @test albs[2] == maximum(albs)  # Maximum albedo for -X, which is facing the Earth

        # Trial 2 (in km, steps)
        rsun = SVector{3, Float64}(1e10   * [0.0, 1.0, 0.0])  
        rsat = SVector{3, Float64}(7e6 * [0.0, 1.0, 0.0])

        data = earth_albedo(rsat, rsun, alb.refl.data)  
        albs = [Simulator.compute_diode_albedo(data, cce, ns[i], rsat) for i = 1:6]
        @test albs[3] == 0.0            # Diode is facing the sun, not Earth
        @test albs[4] == maximum(albs)  # Maximum albedo for -X, which is facing the Earth

        # Trial 3 (new Eam0)
        rsun = SVector{3, Float64}(1e10 * [1.0, 1.0, 1.0])  
        rsat = SVector{3, Float64}(7e6  * [1.0, 1.0, 1.0])

        data = earth_albedo(rsat, rsun, alb.refl.data)  
        albs = [Simulator.compute_diode_albedo(data, cce, ns[i], rsat) for i = 1:6]
        @test all(albs[1:2:end] .== 0.0)
        @test all(albs[2:2:end] .!= 0.0)

        # # Higher at poles than at equator
        rsun = SVector{3, Float64}(1e10 * [0.0, 0.0, 1.0])  
        rsat = SVector{3, Float64}(7e6  * [0.0, 0.0, 1.0])

        data = earth_albedo(rsat, rsun, alb.refl.data)  
        albs_pole = [Simulator.compute_diode_albedo(data, cce, ns[i], rsat) for i = 1:6]

        rsun = SVector{3, Float64}(1e10 * [1.0, 0.0, 0.0])  
        rsat = SVector{3, Float64}(7e6  * [1.0, 0.0, 0.0])

        data = earth_albedo(rsat, rsun, alb.refl.data)  
        albs_equ = [Simulator.compute_diode_albedo(data, cce, ns[i], rsat) for i = 1:6]

        @test norm(albs_pole) > norm(albs_equ)

        # # Higher when in line than when not
        rsun = SVector{3, Float64}(1e10 * [0.0, 1.0, 0.0])  
        rsat = SVector{3, Float64}(7e6  * [0.0, 1.0, 0.0])

        data = earth_albedo(rsat, rsun, alb.refl.data)  
        albs_inline = [Simulator.compute_diode_albedo(data, cce, ns[i], rsat) for i = 1:6]
 
        rsat = SVector{3, Float64}(7e6  * [1.0, 0.0, 0.0])

        data = earth_albedo(rsat, rsun, alb.refl.data)  
        albs_perp = [Simulator.compute_diode_albedo(data, cce, ns[i], rsat) for i = 1:6]

        @test norm(albs_inline) > norm(albs_perp)

    end;

    # Not sure if this is the right gaussian, but it does look gaussian
    @testset "  Rotation Noise" begin 

        @test Simulator.rotation_noise(0.0,  0.0) == I(3)
        @test Simulator.rotation_noise(0.0, 10.0) == I(3)

        @test size(Simulator.rotation_noise(0.1, 1.0)) == (3, 3)

        N = 10000
        v  = [0.0; 1.0; 0.0]
        vs = zeros(N, 3)
        for i = 1:N
            R = Simulator.rotation_noise(0.2, 1.0)
            vs[i, :] .= R * v; 
        end

       
        scatter(vs[:, 1], vs[:, 2], vs[:, 3], title = "Noisy Vectors", 
                    xlim = [-1.5, 1.5], ylim = [-1.5, 1.5], zlim = [-1.5, 1.5],
                    markersize = 1)
        # display( scatter!([0.0], [0.0], [0.0]) )

        mags = [norm(vs[i, :]) for i = 1:N]
        @test all(mags .≈ 1.0)
    end;

    @testset "  Diode Measurement" begin 
        sat = Simulator.SATELLITE()
        alb = get_albedo() 
        x   = Simulator.STATE()
        ecl = 1.0
        sᴵ  = SVector{3, Float64}(1e10 * [1, 0, 0])
        sᴮ  = rand(3); sᴮ = SVector{3, Float64}(sᴮ)

        # Type, Allocs 
        # @code_warntype diode_measurement(sat, alb, x, ecl, sᴵ, sᴮ)

        # Test downsampling
        a₁ = get_albedo(1)
        a₂ = get_albedo(2)
        a₃ = get_albedo(4) 
        sat = Simulator.SATELLITE() 
        x   = Simulator.STATE() 
        ecl = 1.0 
        sᴵ  = SVector{3, Float64}(1e11 * [1.0, 0.0, 0.0])
        sᴮ  = [1.0, 0.0, 1.0];  sᴮ = SVector{3, Float64}(sᴮ / norm(sᴮ))

        function albTest(a₁, a₂, a₃, sat, x, ecl, sᴵ, sᴮ; verbose = false)
            if verbose
                @btime Simulator.diode_measurement($sat, $a₁, $x, $ecl, $sᴵ, $sᴮ);
                @btime Simulator.diode_measurement($sat, $a₂, $x, $ecl, $sᴵ, $sᴮ);
                @btime Simulator.diode_measurement($sat, $a₃, $x, $ecl, $sᴵ, $sᴮ);
            end

            d₁ = Simulator.diode_measurement(sat, a₁, x, ecl, sᴵ, sᴮ);
            d₂ = Simulator.diode_measurement(sat, a₂, x, ecl, sᴵ, sᴮ);
            d₃ = Simulator.diode_measurement(sat, a₃, x, ecl, sᴵ, sᴮ);

            return d₁, d₂, d₃
        end

        d₁, d₂, d₃ = albTest(a₁, a₂, a₃, sat, x, ecl, sᴵ, sᴮ; verbose = false)
        @test d₁[1] ≈ d₂[1] atol = 1e-3

        # Test that I = Ĩ - η, I ≥ 0.0, I = 0.0 during eclipse, etc... 
        I, Ĩ, ηI = Simulator.diode_measurement(sat, a₁, x, ecl, sᴵ, sᴮ);
        @test I[I .> 0.0] ≈ (Ĩ - ηI)[I .> 0.0]   # because we clip at zero, we have to ignore those terms 
        @test all(I .≥ 0.0)

        x = Simulator.STATE()
        sat = Simulator.SATELLITE()
        sᴵ  = SVector{3, Float64}(1e11 * [-0.1, 1.0, 5.0])
        ecl = 0.5
        I, Ĩ, ηI = Simulator.diode_measurement(sat, a₂, x, ecl, sᴵ, sᴮ);
        @test I[I .> 0.0] ≈ (Ĩ - ηI)[I .> 0.0]
        @test all(I .≥ 0.0)

        ecl = 0.0
        I, Ĩ, ηI = Simulator.diode_measurement(sat, a₂, x, ecl, sᴵ, sᴮ);
        @test sum(abs.(I)) == 0.0
        @test Ĩ == ηI 

        # This should also generate an eclipse scenario...
        x = Simulator.STATE(; r = SVector{3, Float64}(7e6 * [-1.0, 0.0, 0.0]))
        sᴵ = SVector{3, Float64}(1e11 * [1.0, 0.0, 0.0])
        ecl = eclipse_cylindrical(vcat([x.r;]...), vcat([sᴵ;]...))
        sᴮ  = x.r / norm(x.r)
        I, Ĩ, ηI = Simulator.diode_measurement(sat, a₂, x, ecl, sᴵ, sᴮ);
        @test sum(abs.(I)) == 0.0
        @test Ĩ == ηI 

        x = Simulator.STATE()
        ecl = 1.0 
        I, Ĩ, ηI = Simulator.diode_measurement(sat, a₂, x, ecl, sᴵ, sᴮ; σ_scale = 0.0);
        @test I == Ĩ 
        @test sum(abs.(ηI)) == 0.0

        # Different number of diodes
        sat = Simulator.SATELLITE(; dio = Simulator.DIODES(; N = 10))
        ecl = 1.0
        I, Ĩ, ηI = Simulator.diode_measurement(sat, a₂, x, ecl, sᴵ, sᴮ; σ_scale = 0.5);
        @test I[I .> 0.0] ≈ (Ĩ - ηI)[I .> 0.0]
        @test all(I .≥ 0.0)

        ecl = 0.0
        I, Ĩ, ηI = Simulator.diode_measurement(sat, a₂, x, ecl, sᴵ, sᴮ);
        @test sum(abs.(I)) == 0.0
        @test Ĩ == ηI 
  
    end;

    @testset "  Sun Measurement" begin 
        
        t = Epoch(2021, 9, 22)  # sun is at ≈[-1.5e11, 1.5e9, 6.9e8]
        dt = 1.0

        # @code_warntype sun_measurement(x, Q, t, dt)

        # In eclipse
        x = Simulator.STATE(; r = SVector{3, Float64}(6.4e6 * [1.0, 0.0, 0.0]))
        Q = SMatrix{3, 3, Float64, 9}(H' * L(x.q) * R(x.q)' * H)
        sᴵ, sᴮ, ecl = Simulator.sun_measurement(x, Q, t)
        @test ecl == 0.0            # In eclipse 
        @test sum(abs.(sᴮ)) != 0.0
        @test sun_position(t) - sᴵ == x.r

        x = Simulator.STATE(; r = SVector{3, Float64}(6.375e6 * [0.0, -1.0, 0.0]))
        Q = SMatrix{3, 3, Float64, 9}(H' * L(x.q) * R(x.q)' * H)
        sᴵ, sᴮ, ecl = Simulator.sun_measurement(x, Q, t)
        @test ecl == 0.0            # In eclipse 
        @test sum(abs.(sᴮ)) != 0.0
        @test sun_position(t) - sᴵ == x.r

        x = Simulator.STATE(; r = SVector{3, Float64}(6.371e3 * [1.0, -1.0, -1.0])) # Bring closer to Earth...
        Q = SMatrix{3, 3, Float64, 9}(H' * L(x.q) * R(x.q)' * H)
        sᴵ, sᴮ, ecl = Simulator.sun_measurement(x, Q, t)
        @test ecl == 0.0            # In eclipse 
        @test sum(abs.(sᴮ)) != 0.0
        

        # In sun
        t = Epoch(2021, 3, 21)  # sun is at ≈[1.5e11, 1.8e9, 7.6e8]

        x = Simulator.STATE(; r = SVector{3, Float64}(6.6e6 * [1.0, 0.0, 0.0]))
        Q = SMatrix{3, 3, Float64, 9}(H' * L(x.q) * R(x.q)' * H)
        sᴵ, sᴮ, ecl = Simulator.sun_measurement(x, Q, t)
        @test ecl == 1.0            # In sun 
        @test sum(abs.(sᴮ)) != 0.0
        @test sun_position(t) - sᴵ == x.r
        @test norm(sᴮ) ≈ 1.0

        x = Simulator.STATE(; r = SVector{3, Float64}(7e6 * [1.0, 1.0, 1.0]))
        Q = SMatrix{3, 3, Float64, 9}(H' * L(x.q) * R(x.q)' * H)
        sᴵ, sᴮ, ecl = Simulator.sun_measurement(x, Q, t)  
        @test ecl == 1.0            # In sun 
        @test sᴮ ≈ Q * (sun_position(t) - x.r) / norm((sun_position(t) - x.r))
        @test norm(sᴮ) ≈ 1.0


        # On corner <- it is very hard to get this 
    end;

    @testset "  Pos Measurement" begin 
        # This testset is pretty boring

        # @code_warntype pos_measurement(x)
        for i = 1:10
            x =  Simulator.STATE();
            r, r̃, ηr = Simulator.pos_measurement(x)
            @test r == x.r 
            @test r ≈ (r̃ - ηr)
        end

        x =  Simulator.STATE();
        r, r̃, ηr = Simulator.pos_measurement(x; σ = 0.0)
        @test r == r̃ 
        @test sum(abs.(ηr)) == 0.0

    end;

    @testset "  Gyro Measurement" begin
        # Also boring 
        
        x = Simulator.STATE()
        # @code_warntype gyro_measurement(x)
        w, w̃, ηw = Simulator.gyro_measurement(x) 
        @test x.ω == w 
        @test w ≈ (w̃ - ηw - x.β) 

        ω = SVector{3, Float64}([1.0, -1.0, 0.1])
        x = Simulator.STATE(ω = ω )
        w, w̃, ηw = Simulator.gyro_measurement(x) 
        @test x.ω == w 
        @test ω == w
        @test w ≈ (w̃ - ηw - x.β) 
        
        x = Simulator.STATE(β = SVector{3, Float64}(zeros(3)))
        w, w̃, ηw = Simulator.gyro_measurement(x) 
        @test x.ω == w 
        @test w ≈ (w̃ - ηw) 

        x = Simulator.STATE() 
        w, w̃, ηw = Simulator.gyro_measurement(x; σ_scale = 0.0) 
        @test x.ω == w 
        @test sum(abs.(ηw)) == 0.0

    end;

    @testset "  Mag Measurement" begin 

        # Perfect Magnetometer (unbiased, no need for calibration) 
        perfect_mag = Simulator.MAGNETOMETER( SVector{3, Float64}(1.0, 1.0, 1.0),   # Scale factor 
                                    SVector{3, Float64}(0.0, 0.0, 0.0),   # Non-orthogonality angles
                                    SVector{3, Float64}(0.0, 0.0, 0.0))   # Bias
         

        # Verify rotation doesn't change vector norm
        sat = Simulator.SATELLITE()
        x   = Simulator.STATE() 
        t   = Epoch(2019, 5, 29)
        dt  = 0.1
        Bᴵ, Bᴮ, B̃ᴮ = Simulator.mag_measurement(sat, x, quat2rot(x.q), t, dt)
        @test norm(Bᴵ) ≈ norm(Bᴮ)

        # Verify a perfect magnetometer results in an unchanged mag vector 
        sat = Simulator.SATELLITE(mag = perfect_mag)
        x   = Simulator.STATE() 
        t   = Epoch(2020, 9, 2)
        dt  = 0.5
        Bᴵ, Bᴮ, B̃ᴮ = Simulator.mag_measurement(sat, x, quat2rot(x.q), t, dt; σ = 0.0)
        @test norm(Bᴵ) ≈ norm(Bᴮ)
        @test Bᴮ == B̃ᴮ

        # Verify that an unrotated sensor has same body/inertial frame 
        sat = Simulator.SATELLITE()
        x   = Simulator.STATE(; q = SVector{4, Float64}(1, 0, 0, 0)) 
        t   = Epoch(2020, 9, 2)
        dt  = 0.75
        Bᴵ, Bᴮ, B̃ᴮ = Simulator.mag_measurement(sat, x, quat2rot(x.q), t, dt)
        @test norm(Bᴵ) ≈ norm(Bᴮ)
        @test Bᴮ != B̃ᴮ
        @test Bᴵ == Bᴮ

    end;

    @testset "  Generate Measurements" begin
        perfect_mag = Simulator.MAGNETOMETER( SVector{3, Float64}(1.0, 1.0, 1.0),   # Scale factor 
                                    SVector{3, Float64}(0.0, 0.0, 0.0),   # Non-orthogonality angles
                                    SVector{3, Float64}(0.0, 0.0, 0.0));  # Bias


        # @code_warntype generate_measurements(sat, alb, x, t, dt)
        # @btime generate_measurements($sat, $alb, $x, $t, $dt);

        # Check with/without noise
        sat = Simulator.SATELLITE();
        alb = get_albedo(2);
        x   = Simulator.STATE(r = SVector{3, Float64}(6.8e6 * [-1.0, 1.0, 0.0]));
        t   = Epoch(2020, 7, 23);  # [-7.7e10, 1e11, 5.2e10]
        dt  = 1.0;
        truth, sensors, ecl, noise = Simulator.generate_measurements(sat, alb, x, t, dt);
        @test ecl == 1.0
        @test truth.t == t 
        @test norm(truth.Bᴵ) ≈ norm(truth.Bᴮ)
        @test all(truth.I .≥ 0)

        # Bias/Noise free system
        sat = Simulator.SATELLITE(; mag = perfect_mag)
        x   = Simulator.STATE(r = SVector{3, Float64}(6.8e6 * [-1.0, 1.0, 0.0]), β = SVector{3, Float64}(zeros(3)) )
        truth, sensor, ecl, noise = Simulator.generate_measurements(sat, alb, x, t, dt;
                                        σB = 0.0, σ_gyro_scale = 0.0, σr = 0.0, σ_current_scale = 0.0);
        @test all(noise.diodes .== 0);
        @test all(noise.gyro .== 0);
        @test all(noise.pos .== 0);
        @test sensor.magnetometer == truth.Bᴮ
        @test sensor.diodes == truth.I
        @test sensor.gyro == x.ω
        @test sensor.pos  == x.r

        #   sun should output ecl = 0.0, so I == 0
        sat = Simulator.SATELLITE();
        alb = get_albedo(2);
        x   = Simulator.STATE(r = SVector{3, Float64}(6.75e6 * [1.0, 0.01, 0.01]));
        t   = Epoch(2019, 11, 3);  # [-1.1e11, -8.8e10, -3.8e10]
        dt  = 5.0;
        truth, sensors, ecl, noise = Simulator.generate_measurements(sat, alb, x, t, dt);
        @test ecl == 0.0
        @test all(truth.I .== 0.0)

        @test truth.t == t 
        @test norm(truth.Bᴵ) ≈ norm(truth.Bᴮ)
        @test all(truth.I .≥ 0)
    end;
end;

