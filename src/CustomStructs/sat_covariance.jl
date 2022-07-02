# [src/MissionSim/CustomStructs/sat_covariance.jl]

# MAY BE UNNECESSARY

""" Originally made so I could define functions on it, but I am not sure it is worht it..."""

struct SAT_COVARIANCE{T}
    """ Kept together b/c may have cross terms """
    """ Probably lower triangular, but not necessarily forced to be so """
    # Too large for Static, (> 100)
    Σ::Matrix{T}   # Covariance matrix for satellite state
    N::Int         # Number of diodes

    function SAT_COVARIANCE(ϕ, β, C, α, ϵ)
        N = size(C, 1)
        Tt = typeof(ϕ[1])

        ℓ = 6 + 3 * N 
        Σ = zeros(Tt, ℓ, ℓ)
        Σ[1:3, 1:3] .= ϕ
        Σ[4:6, 4:6] .= β 

        i₀ = 7; i = i₀ - 1 + N
        Σ[i₀:i, i₀:i]   .= C 

        i₀ = i + 1; i = i₀ - 1 + N
        Σ[i₀:i, i₀:i] .= α 

        i₀ = i + 1; i = i₀ - 1 + N
        Σ[i₀:i, i₀:i] .= ϵ

        Σ = cholesky(Hermitian(Matrix(Σ))).U

        new{Tt}(Σ, N)
    end

    ## Not sure if noise is right!
    function SAT_COVARIANCE(; σϕ = deg2rad(10), σβ = deg2rad(10), σC = 0.1, σα = deg2rad(3), σϵ = deg2rad(3), N = 6)
        """ random """
        Σϕ = diagm( (σϕ^2) * ones(3) )
        Σβ = diagm( (σβ^2) * ones(3) )
        ΣC = diagm( (σC^2) * ones(N) )
        Σα = diagm( (σα^2) * ones(N) )
        Σϵ = diagm( (σϵ^2) * ones(N) )
        
        SAT_COVARIANCE(Σϕ, Σβ, ΣC, Σα, Σϵ)
    end
end
function ϕ(cov::SAT_COVARIANCE{T}) where {T}
    return SMatrix{4, 4, T, 16}(cov.Σ[1:3, 1:3])
end
function β(cov::SAT_COVARIANCE{T}) where {T}
    return SMatrix{3, 3, T, 9}(cov.Σ[4:6, 4:6])
end
function C(cov::SAT_COVARIANCE{T}) where {T}
    N = cov.N
    i₀ = 4
    i  = i₀ + N - 1
    return SMatrix{N, N, T, N * N}(cov.Σ[i₀:i, i₀:i])
end
function α(cov::SAT_COVARIANCE{T}) where {T}
    N = cov.N
    i₀ = 4 + N
    i  = i₀ + N - 1
    return SMatrix{N, N, T, N * N}(cov.Σ[i₀:i, i₀:i])
end
function ϵ(cov::SAT_COVARIANCE{T}) where {T}
    N = cov.N
    i₀ = 4 + 2 * N
    i  = i₀ + N - 1
    return SMatrix{N, N, T, N * N}(cov.Σ[i₀:i, i₀:i])
end
