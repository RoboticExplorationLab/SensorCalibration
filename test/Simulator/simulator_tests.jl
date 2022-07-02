# [src/Simulator/simulator_tests.jl]

""" To Do:
  - to avoid scope conflicts, I am just calling CustomStructs through Simulator rn...
"""

##### Stuff needed to run this independently ####################################################
# using Test, BenchmarkTools
# using StaticArrays, LinearAlgebra, Plots, SatelliteDynamics, EarthAlbedo, JLD2, Distributions

# include("../../src/CustomStructs.jl");  using .CustomStructs 
# include("../../src/Simulator/dynamics.jl")
# include("../../src/Simulator/measurements.jl")
# include("../../src/quaternions.jl")
# include("../../src/mag_field.jl")
#################################################################################################


@testset "Simulator Tests" begin 

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

    #######

    @testset "  Standard Call" begin
        # Call the primary functions
        u = SVector{3, Float64}( zeros(3) ) 
        alb = get_albedo(2);
        sat = SATELLITE();
        x   = STATE();
        t   = Epoch(2021, 12, 25); 
        dt  = 1.0;

        x = rk4(sat.J, x, u, t, dt)
        truth, sensors, ecl, noise = generate_measurements(sat, alb, x, t, dt);
    end;

    @testset "  Sequence, simplified" begin 
        # No noise, bias, fancy gravity terms, etc...
        perfect_mag = MAGNETOMETER(; ideal = true);  

        x = STATE(; β = SVector{3, Float64}(0.0, 0.0, 0.0), ω = SVector{3, Float64}(-0.1, 0.1, 0.25)) 
        u = SVector{3, Float64}( zeros(3) ) 
        alb = get_albedo(2);
        sat = SATELLITE(mag = perfect_mag);
        t   = Epoch(2020, 7, 4); 
        dt  = 2.0;
        N   = 5000;

        truth, sensor, ecl, noise = generate_measurements(sat, alb, x, t, dt);

        truths  = [truth  for i = 1:N];
        sensors = [sensor for i = 1:N];
        ecls    = [ecl    for i = 1:N];
        noises  = [noise  for i = 1:N];
        states  = [x      for i = 1:N];

        for i = 1:N 
            x = rk4(sat.J, x, u, t + (i - 1) * dt, dt; coef_drag = 0.0, coef_srp = 0.0, n_grav = 0, m_grav = 0, third_body = false)
            tr, s, e, n = generate_measurements(sat, alb, x, t + (i - 1) * dt, dt;
                                                σB = 0.0, σ_gyro = 0.0, σr = 0.0, σ_current = 0.0);
            truths[i]  = tr 
            sensors[i] = s 
            ecls[i]    = e 
            noises[i]  = n 
            states[i]  = x
        end

        @test sum([norm(states[i].q) ≈ 1.0  for i = 1:N]) == N
        @test sum([sum(noises[i].diodes) for i = 1:N]) == 0.0

        # plot(truths)      # Mag should be lowest when near equator-ish 
        # plot(sensors)
        # plot(states)
    end;

    @testset "  Sequence, full" begin 

        x = STATE() 
        u = SVector{3, Float64}( zeros(3) ) 
        alb = get_albedo(2);
        sat = SATELLITE();
        t   = Epoch(2020, 7, 4); 
        dt  = 2.0;
        N   = 5000;

        truth, sensor, ecl, noise = generate_measurements(sat, alb, x, t, dt);

        truths  = [truth  for i = 1:N];
        sensors = [sensor for i = 1:N];
        ecls    = [ecl    for i = 1:N];
        noises  = [noise  for i = 1:N];
        states  = [x      for i = 1:N];

        for i = 1:N 
            x = rk4(sat.J, x, u, t + (i - 1) * dt, dt)
            tr, s, e, n = generate_measurements(sat, alb, x, t + (i - 1) * dt, dt);
            truths[i]  = tr 
            sensors[i] = s 
            ecls[i]    = e 
            noises[i]  = n 
            states[i]  = x
        end

        @test sum([norm(states[i].q) ≈ 1.0  for i = 1:N]) == N
        
        # plot(truths)
        # plot(sensors)
        # plot(states)
    end;

end;
