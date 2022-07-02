using Infiltrator, PyCall, Test, MAT
using SatelliteDynamics, LinearAlgebra
using CoordinateTransformations
include("/home/benjj/.julia/dev/EarthAlbedo.jl/src/EarthAlbedo.jl");  using .EarthAlbedo;  __init_albedo__()

function __init__()
    py"""
    import sys
    sys.path.insert(0, "./") # Lets me import local files

    import numpy as np
    from numpy import arccos as acos # update

    from earthAlbedo import EarthAlbedo
    import satelliteDynamics 
    import scipy.io
    

    def sGEODtoECEF(v, flag):
        return satelliteDynamics.sGEODtoECEF(v, use_degrees = flag)

    def spherical2cartesian(v):
        return satelliteDynamics.spherical2cartesian(v)

    def cartesian2spherical(v):
        return satelliteDynamics.cartesian2spherical(v)
    """
end

__init__()


refl_dict = matread("refl.mat") 
model = py"EarthAlbedo"(refl_dict["data"], refl_dict["type"], refl_dict["start_time"], refl_dict["stop_time"])


test_cell_centers = false 
test_gridangle = false
test_cellarea = false
test_earthfov = false
test_idx2rad = false
test_rad2idx = false
test_albedo = true
test_diode = true 

test_sphere2cart = false 
test_cart2sphere = false
test_sGEOD = false

test_cell_centers, test_gridangle, test_cellarea, test_earthfov, test_idx2rad, test_rad2idx, test_albedo, test_diode = true, true, true, true, true, true, true, true 
test_sphere2cart, test_cart2sphere, test_sGEOD = true, true, true


if test_idx2rad
    @testset "idx2rad" begin  # should julias version be exactly one offset?
        sy, sx = 180, 288 

        for i = 1:sy
            for j = 1:sx
                θ, ϵ = idx2rad(i, j, sy, sx)
                theta, eps = model.idx2rad(i-1, j-1, sy, sx) 
                @test ((θ == theta) && (ϵ == eps))
            end
        end
    end
end

if test_rad2idx
    @testset "rad2idx" begin  
        sy, sx = 180, 288 

        for i = 1:sy 
            for j = 1:sx 
                θ, ϵ = idx2rad(i, j, sy, sx)
                î, ĵ = model.rad2idx(θ, ϵ, sy, sx)
                @test ((i == (î+1)) && (j == (ĵ+1)))
            end
        end
    end
end

if test_gridangle
    @testset "gridangle" begin 
        sy, sx = 180, 288 

        i1, i2 = 10, 10
        j1, j2 = 50, 50 
        ρ = gridangle(i1, j1, i2, j2, sy, sx)
        rho = model.gridangle(i1-1, j1-1, i2-1, j2-1, sy, sx)
        @test ((ρ == 0) && (rho == 0))

        for i1 = 1:sy 
            for j1 = 1:sx
                i2, j2 = rand(1:sy), rand(1:sx)
                ρ = gridangle(i1, j1, i2, j2, sy, sx)
                rho = model.gridangle(i1-1, j1-1, i2-1, j2-1, sy, sx)
                @test ρ ≈ rho
            end
        end
    end
end

if test_cellarea
    @testset "cellarea" begin 
        sy, sx = 180, 288 

        for k = 1:100 
            i, j = rand(1:sy), rand(1:sx)

            A = cellarea(i, j, sy, sx)
            A2 = model.cellarea(i-1, j-1, sy, sx)
            @test A ≈ A2
        end  
    end
end

if test_earthfov
    @testset "Earth FOV" begin
        sy, sx = 180, 288  

        # r < Re
        r = 3000   # km 
        θ = deg2rad(rand(1:360))  # radians 
        ϵ = deg2rad(rand(-90:90)) # radians
        fov = earthfov([θ, ϵ, r], sy, sx)
        fov2 = model.earthfov([θ, ϵ, r], sy, sx)
        @test fov ≈ fov2

        # Boundaries
        r = 9000   # km 
        θ = 360    # radians 
        ϵ = deg2rad(rand(-90:90)) # radians
        fov = earthfov([θ, ϵ, r], sy, sx)
        fov2 = model.earthfov([θ, ϵ, r], sy, sx)
        @test fov ≈ fov2

        r = 9000 # km 
        θ = deg2rad(rand(1:360))    # radians 
        ϵ = -90  # radians
        fov = earthfov([θ, ϵ, r], sy, sx)
        fov2 = model.earthfov([θ, ϵ, r], sy, sx)
        @test fov ≈ fov2    

        for k = 1:300
            r = rand(7000:12000)      # km 
            θ = deg2rad(rand(1:360))  # radians 
            ϵ = deg2rad(rand(-90:90)) # radians

            fov = earthfov([θ, ϵ, r], sy, sx)
            fov2 = model.earthfov([θ, ϵ, r], sy, sx)
            @test fov ≈ fov2
        end
    end
end

if test_sGEOD
    @testset "sGEODtoECEF" begin
        for i = 1:100
            lat = rand(-90:90)
            lon = rand(-180:180)
            alt = rand(7000:20000)
            geod = [lon, lat, alt]
            ecef = sGEODtoECEF(geod, use_degrees=true)
            ecef2 = py"sGEODtoECEF"(geod, true)
            @test ecef ≈ ecef2 

            lat = deg2rad(rand(-90:90))
            lon = deg2rad(rand(-180:180))
            alt = rand(7000:20000)
            geod = [lon, lat, alt]
            ecef = sGEODtoECEF(geod, use_degrees = false)
            ecef2 = py"sGEODtoECEF"(geod, false)
            @test ecef ≈ ecef2 
        end
    end
end

if test_cell_centers
    @testset "get_albedo_cell_centers" begin
        cell_centers_ecef = get_albedo_cell_centers()
        cell_centers_alt = model.get_albedo_cell_centers()

        @test cell_centers_ecef ≈ cell_centers_alt
    end
end

if test_diode
    @testset "Diode Albedo" begin
        sy, sx = 180, 288
        albedo_matrix = randn(sy, sx)
        albedo_matrix[abs.(albedo_matrix) .< 0.1] .= 0 

        cell_centers_ecef = get_albedo_cell_centers()

        surf_norm = [1.0, 0.0, 0.0]
        sat = [0.0, 1.0, 0.0]
        diode = get_diode_albedo(albedo_matrix, cell_centers_ecef, surf_norm, sat)
        diode2 = model.get_diode_albedo(albedo_matrix, cell_centers_ecef, surf_norm, sat)
        @test ((diode != 0.0) && (diode2 != 0.0))

        for j = 1:4
            albedo_matrix = randn(sy, sx)
            albedo_matrix[abs.(albedo_matrix) .< 0.1] .= 0 
            for i = 1:50
                surf_norm = randn(3); surf_norm /= norm(surf_norm)
                sat = randn(3); sat /= norm(sat)
                diode  = get_diode_albedo(albedo_matrix, cell_centers_ecef, surf_norm, sat)
                diode2 = model.get_diode_albedo(albedo_matrix, cell_centers_ecef, surf_norm, sat)
                @test diode ≈ diode2
            end
        end    
    end
end

if test_cart2sphere
    Base.isapprox(p1::Spherical, p2) = (p1.r ≈ p2[1]) && (p1.θ ≈ p2[2]) && (p1.ϕ ≈ p2[3])

    @testset "Cart 2 Sphere" begin
        cart = [1.0, 0.0, 0.0]
        sph = py"cartesian2spherical"(cart)
        @test sph == [1.0, 0.0, 0.0]

        for i = 1:100
            cart = [rand(-100:100), rand(-100:100), rand(-100:100)]
            sph = SphericalFromCartesian()(cart)
            sph2 = py"cartesian2spherical"(cart)
            @test sph ≈ sph2
        end

        for i = 1:50 
            cart = [rand(-100:100), rand(-100:100), rand(-100:100)]
            sphere = py"cartesian2spherical"(cart)
            cart_est = py"spherical2cartesian"(sphere)
            @test cart ≈ cart_est
        end
    end
end

if test_sphere2cart
    function correct(vec)
        # Because sometimes its the same angle but ±2π
        while vec[2] < 0.0
            vec[2] += 2*pi
        end
        while vec[2] > 2*pi
            vec[2] -= 2*pi
        end

        return vec
    end

    @testset "Sphere 2 Cart" begin
        sphere = [1.0, deg2rad(90), deg2rad(0)]
        cart = py"spherical2cartesian"(sphere)
        @test cart ≈ [0.0, 1.0, 0.0]

        for i = 1:200
            sphere = [rand(1:100), deg2rad(rand(0:360)), (pi/2) - deg2rad(rand(0:180))]
            spherical = Spherical(sphere[1], sphere[2], sphere[3])
            cart = CartesianFromSpherical()(spherical)
            cart2 = py"spherical2cartesian"(sphere)
            @test cart ≈ cart2
        end

        for i = 1:200
            sphere = [rand(1:100), deg2rad(rand(0:360)), (pi/2) - deg2rad(rand(0:180))]
            cart = py"spherical2cartesian"(sphere)
            sph_est = py"cartesian2spherical"(cart)
            @test correct(sphere) ≈ correct(sph_est)
        end
    end
end

if test_albedo
    @testset "Albedo" begin
        refl_dict = matread("refl.mat")     # Same as ^ but saved in local directory
        sat = [0, 0,  10000000]
        sun = [0, 0, -10000000]
        data = refl_dict["data"]
        refl = refl_struct(refl_dict["data"], refl_dict["type"], refl_dict["start_time"], refl_dict["stop_time"])


        albedos, uni = albedo(sat, sun, refl)
        alb2, un2  = model.albedo(sat, sun, data)
        @test ((sum(abs.(albedos)) == 0) && (alb2 == albedos) && (uni == un2))

        for i = 1:100
            sat = 20000000 * randn(3)
            if norm(sat) < 6500e3
                sat *= 5
            end
            sun = 100000000 * randn(3)
            albedos, uni = albedo(sat, sun, refl)
            alb2,  un2 = model.albedo(sat, sun, data)

            if !((albedos ≈ alb2) && (uni == un2))
                @infiltrate 
            end

            @test (albedos ≈ alb2) && (uni == un2)
        end
    end
end
