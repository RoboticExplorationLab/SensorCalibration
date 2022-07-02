# [src/MissionSim/CustomStructs/ground_truth.jl]

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
function RecipesBase.plot(gt::Vector{GROUND_TRUTH{6, T}}; ds = 1, split = false, kwargs...) where {T}
    N = size(gt, 1)
    ts  = [gt[i].t - gt[1].t for i = 1:ds:N];  # Convert to seconds 
    
    # Split apart and format as matrices for plotting
    Bᴵs = [gt[i].Bᴵ for i = 1:ds:N]; Bᴵs = hcat(Bᴵs...)';
    sᴵs = [gt[i].sᴵ for i = 1:ds:N]; sᴵs = hcat(sᴵs...)';
    ŝᴮs = [gt[i].ŝᴮ for i = 1:ds:N]; ŝᴮs = hcat(ŝᴮs...)';
    Bᴮs = [gt[i].Bᴮ for i = 1:ds:N]; Bᴮs = hcat(Bᴮs...)';
    Is  = [gt[i].I  for i = 1:ds:N]; Is  = hcat(Is... )';

    # Generate plots
    pBᴵ = plot(ts, Bᴵs, title = "Inertial Mag Vec (Bᴵ)", xlabel = "Time (s)", ylabel = "Strength (μT)",  label = ["x" "y" "z"]; kwargs...)
    pBᴮ = plot(ts, Bᴮs, title = "Body Mag Vec (Bᴮ)",     xlabel = "Time (s)", ylabel = "Strength (μT)",  label = ["x" "y" "z"]; kwargs...)
    psᴵ = plot(ts, sᴵs, title = "Inertial Sun Vec (sᴵ)", xlabel = "Time (s)", ylabel = "Distance (m)",   label = ["x" "y" "z"]; kwargs...)
    psᴮ = plot(ts, ŝᴮs, title = "Body Sun Vec (ŝᴮ)",     xlabel = "Time (s)", ylabel = "Relative Scale", label = ["x" "y" "z"]; kwargs...)
    pI  = plot(ts, Is,  title = "Currents (I)",          xlabel = "Time (s)", ylabel = "Current (A)"; kwargs...)
    p   = plot()

    if split   # Return each plot separately
        return pBᴵ, pBᴮ, psᴵ, psᴮ, pI
    else       # Return as one big plot
        return plot(pBᴵ, pBᴮ, psᴵ, psᴮ, pI, p, plot_title = "Ground Truth", layout = (3, 2))
    end
end