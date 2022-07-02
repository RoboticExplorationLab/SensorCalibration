# [src/CustomStructs/ground_truth.jl]

"""
    GROUND_TRUTH{N, T} -> t, Bᴵ, sᴵ, ŝᴮ, Bᴮ, i

      Struct containing the ground truth of various things at a 
    given time stamp. Used to validate the performance of the algorithm, 
    and for aid in debugging. 

    Comes with a custom plotting function.
"""
struct GROUND_TRUTH{N, T}
    t::Epoch                        # Time stamp at which this data holds
    Bᴵ::SVector{3, T}               # Magnetic field vector in inertial frame
    sᴵ::SVector{3, T}               # Sat-Sun vector in inertial frame
    ŝᴮ::SVector{3, T}               # Unit Sat-Sun vector in body frame
    Bᴮ::SVector{3, T}               # Magnetic field vector in body frame
    I::SVector{N, T}                # Diode currents 
end

""" plot(gt, s; start, stop, ds, kwargs)

      Plotting recipe for a vector of GROUND_TRUTH structs. Can plot the five data sets over time (which 
    is the default), or just one specific portion, which is selected by the second argument. Keyword arguments 
    allow for using only a portion of the history rather than all. Also includes a downsample rate `ds` 
    for long runs. 

    Available symbols:
      `:a`: Show all sensor measurements
      `:Bᴵ`: Show mag field in inertial frame
      `:sᴵ`: Show sun vector in inertial frame  
      `:sᴮ`: Show sun vector in body frame (unit)
      `:Bᴮ`: Show mag field in body frame
      `:I`:  Show generated current
"""
function RecipesBase.plot(gt::Vector{GROUND_TRUTH{6, T}}, t::Symbol = :a; start = 1, stop = nothing, ds = 2, split = false, kwargs...) where {T}
    N = isnothing(stop) ? size(gt, 1) : stop 


    ts  = [gt[i].t - gt[1].t for i = start:ds:N];  # Convert to seconds 

    if t == :a    
        # Split apart and format as matrices for plotting
        Bᴵs = [gt[i].Bᴵ for i = start:ds:N]; Bᴵs = reduce(hcat, Bᴵs)';
        sᴵs = [gt[i].sᴵ for i = start:ds:N]; sᴵs = reduce(hcat, sᴵs)';
        ŝᴮs = [gt[i].ŝᴮ for i = start:ds:N]; ŝᴮs = reduce(hcat, ŝᴮs)';
        Bᴮs = [gt[i].Bᴮ for i = start:ds:N]; Bᴮs = reduce(hcat, Bᴮs)';
        Is  = [gt[i].I  for i = start:ds:N]; Is  = reduce(hcat, Is)';

        pBᴵ = plot(ts, Bᴵs, title = "Bᴵ", xlabel = "Time (s)", ylabel = "Strength (μT)",  label = ["x" "y" "z"]; kwargs...)
        pBᴮ = plot(ts, Bᴮs, title = "Bᴮ", xlabel = "Time (s)", ylabel = "Strength (μT)",  label = ["x" "y" "z"]; kwargs...)
        psᴵ = plot(ts, sᴵs, title = "sᴵ", xlabel = "Time (s)", ylabel = "Distance (m)",   label = ["x" "y" "z"]; kwargs...)
        psᴮ = plot(ts, ŝᴮs, title = "ŝᴮ", xlabel = "Time (s)", ylabel = "Relative Scale", label = ["x" "y" "z"]; kwargs...)
        pI  = plot(ts, Is,  title = "I",  xlabel = "Time (s)", ylabel = "Current (A)"; kwargs...)
        
        return plot(pBᴵ, pBᴮ, psᴵ, psᴮ, pI, plot_title = "Ground Truth")

    elseif t == :Bᴵ
        Bᴵs = [gt[i].Bᴵ for i = start:ds:N]; Bᴵs = reduce(hcat, Bᴵs)';
        pBs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            pB = plot(ts, Bᴵs[:, 1], xlabel = "Time (s)", ylabel = "Strength (μT)", label = labels[i]; kwargs...)
            push!(pBs, pB)
        end
        return plot(pBs..., plot_title = "Ground Truth - Bᴵ")

    elseif t == :sᴵ 
        sᴵs = [gt[i].sᴵ for i = start:ds:N]; sᴵs = reduce(hcat, sᴵs)';
        pss = []
        labels = ["x" "y" "z"]
        for i = 1:3
            ps = plot(ts, sᴵs[:, i], xlabel = "Time (s)", ylabel = "Distance (m)", label = labels[i]; kwargs...)
            push!(pss, ps)
        end
        return plot(pss..., plot_title = "Ground Truth - sᴵ")

    elseif t == :sᴮ 
        ŝᴮs = [gt[i].ŝᴮ for i = start:ds:N]; ŝᴮs = reduce(hcat, ŝᴮs)';
        pss = []
        labels = ["x" "y" "z"]
        for i = 1:3 
            ps =  plot(ts, ŝᴮs[:, i], xlabel = "Time (s)", ylabel = "Relative Scale", label = labels[i]; kwargs...)
            push!(pss, ps)
        end
        return plot(pss..., plot_title = "Ground Truth - sᴮ (unit)")

    elseif t == :Bᴮ 
        Bᴮs = [gt[i].Bᴮ for i = start:ds:N]; Bᴮs = reduce(hcat, Bᴮs)';
        pbs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            pb = plot(ts, Bᴮs[:, i], title = "Bᴮ", xlabel = "Time (s)", ylabel = "Strength (μT)", label = labels[i]; kwargs...)
            push!(pbs, pb)
        end
        return plot(pbs..., plot_title = "Ground Truth - Bᴮ")

    elseif t == :I 
        Is  = [gt[i].I  for i = start:ds:N]; Is  = reduce(hcat, Is)';
        pis = []
        for i = 1:6 
            pi = plot(ts, Is[:, i],  title = "I", xlabel = "Time (s)", ylabel = "Current (A)"; kwargs...)
            push!(pis, pi)
        end
        return plot(pis..., plot_title = "Ground Truth - I")

    else 
        println("Warning: Type $t is not ready for plotting!")
        println("\tViable symbols for GROUND_TRUTH are: Bᴵ, sᴵ, Bᴮ, sᴮ, I")
    end
end

