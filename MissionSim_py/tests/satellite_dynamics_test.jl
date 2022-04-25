using Infiltrator, Test, PyCall 
using LinearAlgebra, SatelliteDynamics
include("mag_field.jl")


function __init_py__()
    py"""
    import sys 
    sys.path.insert(0, "./")

    import numpy as np 
    from satelliteDynamics import *

    """
end

__init_py__()


if false
    @testset "Epoch" begin

        epc_jl = Epoch(2021, 9, 1, 23, 0, 0, 0.0) 
        epc_py = py"Epoch"(2021, 9, 1, 23, 0, 0, 0.0) 

        @test epc_jl.days == epc_py.days 
        @test epc_jl.seconds == epc_py.seconds 
        @test epc_jl.nanoseconds == epc_py.nanoseconds 
        @test epc_jl.tsys == epc_py.tsys
    end

    @testset "Sun Position" begin
        
        for i = 1:100
            # @show i 
            year = 2021 + rand(-5:10)
            month = rand(1:12)
            day = rand(1:28)
            hour = rand(1:23)
            minute = rand(1:59)
            second = rand(1:59)
            nanosecond = abs(100.0 * randn())

            epc_jl = Epoch(year, month, day, hour, minute, second, nanosecond)
            epc_py = py"Epoch"(year, month, day, hour, minute, second, nanosecond)

            p1 = sun_position(epc_jl)
            p2 = py"sun_position"(epc_py)

            @test p1 ≈ p2
        end 
    end
end


##### TIME stuff 

function get_random_epoch()
    year = rand(2018:2021)
    month = rand(1:12)
    if month == 2
        day = rand(1:28)
    elseif ((month == 4) || (month == 6) || (month == 9) || (month == 11))
        day = rand(1:30)
    else
        day = rand(1:31)
    end
    hour = rand(1:23)
    min = rand(1:59)
    sec = rand(1:59)
    nano = rand(1:(1e9 - 1))

    epc_jl = Epoch(year, month, day, hour, min, sec, nano)
    epc_py = py"Epoch.ymd"(year, month, day, hour, min, sec, nano)

    return epc_jl, epc_py
end

if false
    @testset "Cal Date" begin
        for i = 1:500
            epc_jl, epc_py = get_random_epoch()
            
            y1, m1, d1, h11, h21, h31, h41 = SatelliteDynamics.caldate(epc_jl)
            y2, m2, d2, h12, h22, h32, h42 = py"caldate"(epc_jl)
            @test (((y1 == y2)  && (m1 == m2) && (d1 == d2) && 
                    (h11 == h12.__float__()) && (h21 == h22.__float__()) && 
                    (h31 == h32.__float__()) && (h41 == h42.__float__())))
        end              
    end

    @testset "Day of Year" begin
        for i = 1:500
            epc_jl, epc_py = get_random_epoch()
            
            doy1 = SatelliteDynamics.day_of_year(epc_jl)
            doy2 = py"day_of_year"(epc_jl)
            @test (doy1 == doy2)
        end
    end

    @testset "Caldate 2 MJD" begin

        for i = 1:500
            year, month, day = rand(2015:2022), rand(1:12), rand(1:28)

            mjd = SatelliteDynamics.caldate_to_mjd(year, month, day)
            m2 = py"caldate_to_mjd"(year, month, day)

            @test mjd == m2
        end
    end

    @testset "Modified Julian Date (mjd)" begin
        epc = get_random_epoch();
        mjd_jl = mjd(epc, tsys = "TT")
        mjd_py = py"mjd"(epc, tsys = "TT")
        @test mjd_jl ≈ mjd_py
    end
end



##### IGRF Stuff (not actually SatelliteDynamics)
# TODO add in test for IGRF13, not just my_igrf_13

if false 
    @testset "My IGRF" begin
        for i = 1:500
            date = rand(-5000:5000)
            alt = rand(-100:10000) * 1.0
            lat = rand(-90:90)
            elon = rand(-180:180)
            order = 13

            x,  y,  z  = my_igrf_13(date, alt, lat, elon, order)
            x2, y2, z2 = py"my_igrf_13"(date, alt, lat, elon, order)
            @test ((x ≈ x2) && (y ≈ y2) && (z ≈ z2))
        end
    end

    @testset "ECEF Q NED" begin
        for i = 1:1000
            lat = rand(-90:90)
            elon = rand(-180:180)

            M = ecef_Q_ned_mat(elon, lat)
            M2 = py"ecef_Q_ned_mat"(elon, lat)
            @test (M ≈ M2)
        end
    end
end

##### OrbitDynamics.jl Portion 

# ASSUMES SatelliteDynamics.eclipse_conical HAS NOT BEEN FIXED!!
if false
    @testset "Eclipse Conical" begin 

        pos = [-1e7, 0, 0]
        sun = [7e10,  0, 0]
        ν_jl = eclipse_conical(-pos, sun)
        ν_py = py"eclipse_conical"(-pos, sun)

        @test ((ν_jl == ν_py) && (ν_jl == 0.0))

        pos = -pos
        ν_jl = eclipse_conical(-pos, sun)
        ν_py = py"eclipse_conical"(-pos, sun)

        @test ((ν_jl == ν_py) && (ν_jl == 1.0))

        # Parital 
        pos = [2.0007668461654785e8, -9.345982592623848e6, -1.2521231216893423e8]
        sun = [-5.131017076755866e10, 3.3619960178098483e9, 3.3922438206486088e10]
        ν_jl = eclipse_conical(-pos, sun)
        ν_py = py"eclipse_conical"(-pos, sun)

        @test ((ν_jl ≈ ν_py)) # ∈ (0, 1)


        for i = 1:100000
            pos, sun = 0, 0
            while (norm(pos) < 7e6) || (norm(sun) < 1e9)
                pos = 1e8 * randn(3)
                sun = 5e10 * randn(3)
            end

            ν_jl = eclipse_conical(-pos, sun)
            ν_py = py"eclipse_conical"(-pos, sun)

            @test (ν_jl ≈ ν_py)
        end
    end
end





