

# Test for duplication === 



@testset "Diode angle wrapping in SAT_STATE" begin

    Nₜ = 1000
    ts = zeros(Nₜ)
    for i = 1:Nₜ
        α = SVector{3, Float64}(5 * randn(3))
        b = CustomStructs.wrap(α) #, -pi, pi)

        ts[i] = all(b .> -pi) && (all(b .≤ pi)) && all(sin.(α) .≈ sin.(b)) && all(cos.(α) .≈ cos.(b))
    end
    @test all(ts .== 1)

    # Test that α and ϵ are unchanged
    ts = zeros(Nₜ)
    for i = 1:Nₜ
        α = deg2rad.(rand(-179:1:180, 6)) 
        ϵ = deg2rad.(rand(-89:1:89,   6))
         
        α₁, ϵ₁ = CustomStructs.cart2sph(CustomStructs.sph2cart(α, ϵ)...)

        ts[i] = all(α .> -pi) && all(α .≤ pi) && all(ϵ .≤ pi/2) && all(ϵ .> -pi/2) && (α ≈ α₁) && (ϵ ≈ ϵ₁)
    end
    @test all(ts .== 1)


    # Test that α and ϵ result in the same Cartesian position 
    ts = zeros(Nₜ)
    for i = 1:Nₜ
        α = 5 * randn(6)
        ϵ = 5 * randn(6)
        x₀, y₀, z₀ = CustomStructs.sph2cart(α, ϵ)

        α₁, ϵ₁ = CustomStructs.cart2sph(x₀, y₀, z₀)
        x₁, y₁, z₁ = CustomStructs.sph2cart(α₁, ϵ₁)

        ts[i] = all(α₁ .> -pi) && all(α₁ .≤ pi) && all(ϵ₁ .≤ pi/2) && all(ϵ₁ .> -pi/2) && all(x₀ ≈ x₁) && all(y₀ ≈ y₁) && all(z₀ ≈ z₁)
    end
    @test all(ts .== 1)

end