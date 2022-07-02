# [src/CustomStructs/noise.jl]

"""
    NOISE{T} -> diodes, gyro, pos 

      Struct containing the noise in each sensor measurement, except for 
    some resaon I don't have magnetometer noise. This is really only used 
    for debugging anyway, so...

"""
struct NOISE{N, T}
    diodes::SVector{N, T}       # Noise in each measured diode current
    gyro::SVector{3, T}         # Noise in the gyroscope measurement
    pos::SVector{3, T}          # Noise in the position estimate 
end

""" plot(n, t; start, stop, ds, kwargs)

      Plotting recipe for a vector of NOISE structs. Can plot all the noise for the three tracked sensor measurements
    over time (which is the default), or just one specific sensor, which is selected by the second argument. Keyword arguments 
    allow for using only a portion of the history rather than all. Additionally, there is a parameter for downsampling.

    Available symbols:
      `:a`: Show all sensor noise
      `:d`: Show magnetometer noise
      `:g`: Show gyroscope noise
      `:p`: Show position noise
"""
function RecipesBase.plot(n::Vector{NOISE{6, T}}, t::Symbol = :a; start = 1, stop = nothing, ds = 1, kwargs...) where {T}
    N = isnothing(stop) ? size(n, 1) : stop
    
    if t == :a 
        dis = [n[i].diodes for i = start:ds:N]; dis = reduce(hcat, dis)';
        gys = [n[i].gyro   for i = start:ds:N]; gys = reduce(hcat, gys)';      
        pos = [n[i].pos    for i = start:ds:N]; pos = reduce(hcat, pos)';

        pd  = plot(dis, title = "Diode Noise",    xlabel = "Index (ds: $ds)", ylabel = "Value", label = ["1" "2" "3" "4" "5" "6"]; kwargs...)
        pg  = plot(rad2deg.(gys), title = "Gyro Noise",     xlabel = "Index (ds: $ds)", ylabel = "Value (deg/s)", label = ["x" "y" "z"]; kwargs...)
        pp  = plot(pos, title = "Position Noise", xlabel = "Index (ds: $ds)", ylabel = "Value", label = ["x" "y" "z"]; kwargs...)
        
        return plot(pd, pg, pp, plot_title = "Noise")

    elseif t == :d 
        dis = [n[i].diodes for i = start:ds:N]; dis = reduce(hcat, dis)';
        pds = []
        labels = ["1" "2" "3" "4" "5" "6"]
        for i = 1:6 
            pd = plot(dis[:, i], xlabel = "Index (ds: $ds)", ylabel = "Value", label = labels[i]; kwargs...)
            push!(pds, pd)
        end
        return plot(pds..., plot_title = "Noise - Diodes")

    elseif t == :g 
        gys = [n[i].gyro   for i = start:ds:N]; gys = reduce(hcat, gys)';   
        pgs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            pg = plot(rad2deg.(gys[:, i]), xlabel = "Index (ds: $ds)", ylabel = "Value (deg/s)", label = labels[i]; kwargs...)   
            push!(pgs, pg)
        end 
        return plot(pgs..., plot_title = "Noise - Gyro")
  
    elseif t == :p 
        pos = [n[i].pos    for i = start:ds:N]; pos = reduce(hcat, pos)';
        pps = []
        labels = ["x" "y" "z"]
        for i = 1:3
            pp = plot(pos[:, i], xlabel = "Index (ds: $ds)", ylabel = "Value", label = labels[i]; kwargs...)
            push!(pps, pp)
        end
        return plot(pps..., plot_title = "Noise - Position")

    else
      println("Warning: Type $t is not ready for plotting!")
      println("\tViable symbols for NOISE are: a, d, g, p")
    end
end

