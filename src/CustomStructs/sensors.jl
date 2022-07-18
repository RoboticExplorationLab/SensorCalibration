# [src/CustomStructs/sensors.jl]

"""
    SENSORS{N, T} -> magnetometer, diodes, gyro, pos

      Struct containing the measurements for each of the satellite sensors 
    (magnetometer, diodes, gyroscope, and something for position). Also 
    comes with a custom plotting recipe. 
"""
struct SENSORS{N, T}
    magnetometer::SVector{3, T}      # Measured magnetometer reading
    diodes::SVector{N, T}            # Measured photodiode currents (per diode)
    gyro::SVector{3, T}              # Measured angular velocity
    pos::SVector{3, T}               # Measured position 
end

""" plot(s, sensor; start, stop, kwargs)

      Plotting recipe for a vector of SENSOR structs. Can plot all four sensor measurements over time (which 
    is the default), or just one specific sensor, which is selected by the second argument. Keyword arguments 
    allow for using only a portion of the history rather than all.

    Available symbols:
      `:a`: Show all sensor measurements
      `:m`: Show magnetometer measurements
      `:d`: Show diode measurements 
      `:g`: Show gyroscope measurements
      `:p`: Show position measurements
"""
function RecipesBase.plot(s::Vector{SENSORS{6, T}}, sensor::Symbol = :a; start = 1, stop = nothing, kwargs...) where {T}
    N = isnothing(stop) ? size(s, 1) : stop

    # Split apart and format as matrices for plotting
    if sensor == :a
        mags = [s[i].magnetometer for i = start:N]; mags = reduce(hcat, mags)'; 
        dios = [s[i].diodes       for i = start:N]; dios = reduce(hcat, dios)'; 
        gyrs = [s[i].gyro         for i = start:N]; gyrs = reduce(hcat, gyrs)'; 
        poss = [s[i].pos          for i = start:N]; poss = reduce(hcat, poss)'; 

        # Make the plots
        pM = plot(mags, title = "Magnetometers", xlabel = "Index", ylabel = "Mag Field (μT)",   label = ["x" "y" "z"]; kwargs...)
        pD = plot(dios, title = "Diodes",        xlabel = "Index", ylabel = "Current (A?)",                         ; kwargs...)
        pG = plot(rad2deg.(gyrs), title = "Gyroscope",     xlabel = "Index", ylabel = "Ang Vel (rad/s)", label = ["x" "y" "z"]; kwargs...)
        pP = plot(poss, title = "Position",      xlabel = "Index", ylabel = "Position (m)",    label = ["x" "y" "z"]; kwargs...)
        return plot(pM, pD, pG, pP, plot_title = "Sensors", layout = (2, 2))
    
    elseif sensor == :m 
        mags = [s[i].magnetometer for i = start:N]; mags = reduce(hcat, mags)';
        ms = []
        for i = 1:3 
            m = plot(mags[:, i])
            push!(ms, m)
        end

        return plot(ms..., plot_title = "Magnetometer Measurements")

    elseif sensor == :d
        dios = [s[i].diodes for i = start:N]; dios = reduce(hcat, dios)'; 
        ds = []
        nd = size(s[1].diodes, 1)
        labels = ["1", "2", "3", "4", "5", "6"]
        for i = 1:nd
            d = plot(dios[:, i], label = labels[i])
            push!(ds, d)
        end
        return plot(ds..., plot_title = "Diode Measurements") 
        
    elseif sensor == :g 
        gyrs = [s[i].gyro for i = start:N]; gyrs = reduce(hcat, gyrs)';
        gs = []
        labels = ["x", "y", "z"]
        for i = 1:3 
            g = plot(rad2deg.(gyrs[:, i]), label = labels[i], xlabel = "Index", ylabel = "Vel (deg/s)")
            push!(gs, g)
        end
        pn = plot(rad2deg.(ng), label = "Magnitude")
        push!(gs, pn)

        return plot(gs..., plot_title = "Gyro Measurements")

    elseif sensor == :p
        poss = [s[i].pos          for i = start:N]; poss = reduce(hcat, poss)';
        ps = []
        labels = ["x", "y", "z"]
        for i = 1:3 
            p = plot(poss[:, i], label = labels[i])
            push!(ps, p)
        end
        return plot(ps..., plot_title = "Position Measurements")

    else 
        println("Warning! Symbol $sensor not supported yet!")
        println("\tAvailable symbols: a, m, d, g, p")
    end
end

