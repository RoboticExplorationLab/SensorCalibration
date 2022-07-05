# [src/MissionSim/CustomStructs/sat_covariance.jl]

# MAY BE UNNECESSARY

""" Originally made so I could define functions on it, but I am not sure it is worht it..."""

struct SAT_COVARIANCE{T}
    """ Kept together b/c may have cross terms """
    """ Probably lower triangular, but not necessarily forced to be so """
    # Too large for Static, (> 100)
    Σ::Matrix{T}   # Covariance matrix for satellite state
    N::Int         # Number of diodes

    function SAT_COVARIANCE(ϕ, βᵧ, s, ζ, βₘ)
        N = 3 #size(s, 1)
        Tt = typeof(ϕ[1])

        ℓ =15
        Σ = zeros(Tt, ℓ, ℓ)
        Σ[1:3, 1:3] .= ϕ
        Σ[4:6, 4:6] .= βᵧ
        Σ[7:9, 7:9] .= s 
        Σ[10:12, 10:12] .= ζ 
        Σ[13:15, 13:15] .= βₘ

        Σ = cholesky(Hermitian(Matrix(Σ))).U

        new{Tt}(Σ, N)
    end

    ## Not sure if noise is right!
    function SAT_COVARIANCE(; σϕ = deg2rad(10), σβᵧ = deg2rad(10), σs = 0.25, σζ = deg2rad(5), σβₘ = deg2rad(6.0))  # 0.1, 3, 4
        """ random """
        Σϕ  = diagm( (σϕ^2) * ones(3) )
        Σβᵧ = diagm( (σβᵧ^2) * ones(3) )
        Σs  = diagm( (σs^2) * ones(3) )
        Σζ  = diagm( (σζ^2) * ones(3) )
        Σβₘ = diagm( (σβₘ^2) * ones(3) )
        
        SAT_COVARIANCE(Σϕ, Σβᵧ, Σs, Σζ, Σβₘ)
    end
end
function ϕ(cov::SAT_COVARIANCE{T}) where {T}
    return SMatrix{4, 4, T, 16}(cov.Σ[1:3, 1:3])
end
function β(cov::SAT_COVARIANCE{T}) where {T}
    return SMatrix{3, 3, T, 9}(cov.Σ[4:6, 4:6])
end