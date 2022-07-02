# [src/CustomStructs/state.jl]

"""
    STATE{S, T} -> r, v, q, ω, β, x 

      Current state of the environment, used for propagating dynamics 
    and generating measurements. Not to be confused with the SATELLITE state
    SAT_STATE. 

    Comes with a custom plotting recipe, as well as a few basic functions 
    for addition, multiplication, and division of states.
"""
struct STATE{T}
    r::SVector{3, T}   # Position (Cartesian)
    v::SVector{3, T}   # Velocity
    q::SVector{4, T}   # Scalar-first quaternion
    ω::SVector{3, T}   # Angular velocity
    β::SVector{3, T}   # Gyroscope bias


    function STATE(r::SVector{3, T}, v::SVector{3, T}, q::SVector{4, T}, ω::SVector{3, T}, β::SVector{3, T}) where {T} 
        """ Primary Constructor """

        # q = (norm(q) == 1) ? q  : (q / norm(q))  # Enforce unit norm constraint
        # q = (q[1] < 0) ? -q : q
        new{T}(r, v, q, ω, β)
    end

    function STATE(r::Vector{T}, v::Vector{T}, q::Vector{T}, ω::Vector{T}, β::Vector{T}) where {T}
        # @warn "Initializing STATE without static arrays! Automatically converting..."
        _r = SVector{3, T}(r)
        _v = SVector{3, T}(v)
        _q = SVector{4, T}(q)
        _ω = SVector{3, T}(ω)
        _β = SVector{3, T}(β)

        STATE(_r, _v, _q, _ω, _β)  # Call the default constructor
    end

    function STATE(; r = nothing, v = nothing, q = nothing, ω = nothing, β = nothing, a = nothing, Rₑ = 6378.1363e3, μ = 3.9860044188e14)
        """ Random Constructor - Primarily used for testing and convenience """
        
        function rand_quat()
            q = randn(4)
            return q / norm(q)
        end

        _r = isnothing(r) ? (Rₑ + (rand(20:100) * 10^4)) * [1.0, 0.0, 0.0] : r  
        _a = isnothing(a) ? norm(_r) * (1.0 + 0.5 * rand()) : a
        _v = isnothing(v) ? sqrt(μ*( (2/norm(_r)) - (1/_a))) * [0.0, 1.0, 0.0] : v 
        _q = isnothing(q) ? rand_quat() : q 
        _ω = isnothing(ω) ? randn(3) : ω
        _β = isnothing(β) ? randn(3) : β

        STATE(SVector{3, Float64}(_r), SVector{3, Float64}(_v), SVector{4, Float64}(_q), 
                SVector{3, Float64}(_ω), SVector{3, Float64}(_β))
    end
end

""" Converts a struct to a vector """
function x(s::STATE{T}) where {T}
    vec = [s.r; s.v; s.q; s.ω; sβ]
    return SVector{16, T}(vec)
end

""" plot(s, sensor; start, stop, kwargs)

      Plotting recipe for a vector of SENSOR structs. Can plot all four sensor measurements over time (which 
    is the default), or just one specific sensor, which is selected by the second argument. Keyword arguments 
    allow for using only a portion of the history rather than all.

    Available symbols:
      `:a`: Show all sensor measurements
      `:r`: Show position
      `:v`: Show velocity
      `:q`: Show attitude (scalar-first quaternion)
      `:ω`: Show angular velocity 
      `:β`: Show gyroscope bias
"""
function RecipesBase.plot(s::Vector{STATE{T}}, t::Symbol = :a; start = 1, stop = nothing, ds = 1, kwargs...) where {T}
    N = isnothing(stop) ? size(s, 1) : stop 

    if t == :a
        # Split apart and format as matrices for plotting
        rs = [vcat([s[i].r;]...) for i = start:ds:N]; rs = hcat(rs...)';
        vs = [vcat([s[i].v;]...) for i = start:ds:N]; vs = hcat(vs...)';
        qs = [vcat([s[i].q;]...) for i = start:ds:N]; qs = hcat(qs...)';  
        ωs = [vcat([s[i].ω;]...) for i = start:ds:N]; ωs = hcat(ωs...)';
        mag_ω = [norm(s[i].ω) for i = start:ds:N]
        βs = [vcat([s[i].β;]...) for i = start:ds:N]; βs = hcat(βs...)';

        pr = plot(rs, title = "Position (Cart)", xlabel = "Index", ylabel = "Position (m)",    label = ["x" "y" "z"]; kwargs...)
        pv = plot(vs, title = "Velocity",        xlabel = "Index", ylabel = "Velocity (m/s)",  label = ["x" "y" "z"]; kwargs...)
        pq = plot(qs, title = "Attitude (quat)", xlabel = "Index",                             label = ["x" "y" "z"]; kwargs...)
        pω = plot(rad2deg.(ωs), title = "Ang Velocity",    xlabel = "Index", ylabel = "Ang Vel (rad/s)", label = ["x" "y" "z"]; kwargs...);
            pω = plot!(mag_ω, label = "||ω||", c = :black, ls = :dash)
        pβ = plot(βs, title = "Gyro Bias",       xlabel = "Index", ylabel = "Bias (rad/s)",    label = ["x" "y" "z"]; kwargs...)
        p   = plot()

        return plot(pr, pv, pq, pω, pβ, p, plot_title = "State", layout = (3, 2))
    elseif t == :r
        pos = [vcat([s[i].r;]...) for i = start:ds:N]; pos = hcat(pos...)';
        rs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            r = plot(pos[:, i], label = labels[i])
            push!(rs, r)
        end
        return plot(rs..., layout = 3, plot_title = "Position (Cart)", xlabel = "Index", ylabel = "Position (m)"; kwargs...)

    elseif t == :v  
        vels = [vcat([s[i].v;]...) for i = start:ds:N]; vels = hcat(vels...)';
        vs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            v = plot(vels[:, i], label = labels[i])
            push!(vs, v)
        end
        return plot(vs..., plot_title = "Velocity", xlabel = "Index", ylabel = "Velocity (m/s)"; kwargs...)

    elseif t == :q 
        quats = [vcat([s[i].q;]...) for i = start:ds:N]; quats = hcat(quats...)';
        qs = []
        for i = 1:4
            q = plot(quats[:, i])
            push!(qs, q)
        end
        return plot(qs..., plot_title = "Attitude (quat)", xlabel = "Index", ylabel = "Magnitutde"; kwargs...)


    elseif t == :ω 
        angs = [vcat([s[i].ω;]...) for i = start:ds:N]; angs = hcat(angs...)';
        ωs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            ω = plot(rad2deg.(angs[:, i]), label = labels[i])
            push!(ωs, ω)
        end
        nω = [norm(rad2deg.(angs[i, :])) for i = start:ds:N]
        push!(ωs, plot(nω, title = "Magnitude"))
        return plot(ωs..., plot_title = "Ang Vel", xlabel = "Index", ylabel = "Vel (deg/s)"; kwargs...)


    elseif t == :β
        bias = [vcat([s[i].β;]...) for i = start:ds:N]; bias = hcat(bias...)';
        βs = []
        labels = ["x" "y" "z"]
        for i = 1:3
            β = plot(bias[:, i], label = "labels[i]")
            push!(βs, β)
        end
        return plot(βs..., plot_title = "Gyro Bias", xlabel = "Index", ylabel = "Magnitutde"; kwargs...)


    else 
        println("Warning: State $t is not ready for plotting!")
        println("\tViable symbols: r, v, q, ω, β, a")
    end

end

""" Normalizes a quaternion """
function normalize_quat(x::STATE)
    q = x.q / norm(x.q)
    return STATE(x.r, x.v, q, x.ω, x.β)
end

# Define addition, multiplication, and subtraction?  (NEEDED for RK4)
function Base.:+(x₁::STATE, x₂::STATE) 
    r = x₁.r + x₂.r 
    v = x₁.v + x₂.v 
    q = x₁.q + x₂.q 
    ω = x₁.ω + x₂.ω 
    β = x₁.β + x₂.β 

    # q /= norm(q)
    return STATE(r, v, q, ω, β)
end

function Base.:*(k::Real, x::STATE)
    r = k * x.r 
    v = k * x.v  
    q = k * x.q  
    ω = k * x.ω  
    β = k * x.β  

    # q /= norm(q)
    return STATE(r, v, q, ω, β)
end

function Base.:/(x::STATE, k::Real)
    r = x.r / k
    v = x.v / k
    q = x.q / k  
    ω = x.ω / k
    β = x.β / k 

    return STATE(r, v, q, ω, β)
end