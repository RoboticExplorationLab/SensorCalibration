# [src/reports.jl]

"""
    Contains functions for evaluating the performance of the Sim/Filter setup 
"""

function magnetometer_calibration_report(sat_true::SATELLITE, est_hist::Vector{SATELLITE{6, T}}) where {T}
    N = size(est_hist, 1)
    s_est, ζ_est, β_est = zeros(3, N), zeros(3, N), zeros(3, N);

    for i = 1:N
        s_est[:, i] = est_hist[i].magnetometer.scale_factors
        ζ_est[:, i] = rad2deg.(est_hist[i].magnetometer.non_ortho_angles)
        β_est[:, i] = est_hist[i].magnetometer.bias
    end

    aᶠ, bᶠ, cᶠ = round.(est_hist[N].magnetometer.scale_factors, sigdigits = 3)
    ρᶠ, λᶠ, νᶠ = round.(rad2deg.(est_hist[N].magnetometer.non_ortho_angles), sigdigits = 3)
    βxᶠ, βyᶠ, βzᶠ = round.(est_hist[N].magnetometer.bias, sigdigits = 3)

    a, b, c = round.(sat_true.magnetometer.scale_factors, sigdigits = 3)
    ρ, λ, ν = round.(rad2deg.(sat_true.magnetometer.non_ortho_angles), sigdigits = 3)
    βx, βy, βz = round.(sat_true.magnetometer.bias, sigdigits = 3)

    a₀, b₀, c₀ = round.(est_hist[1].magnetometer.scale_factors, sigdigits = 3)
    ρ₀, λ₀, ν₀ = round.(rad2deg.(est_hist[1].magnetometer.non_ortho_angles), sigdigits = 3)
    βx₀, βy₀, βz₀ = round.(est_hist[1].magnetometer.bias, sigdigits = 3)

    println("__________________________________________________________________________")
    println("___PARAM___|___Truth____|__Final Guess__|__Init Guess__|__Improved?__")
    println("     a     |   $a\t|    $aᶠ  \t|    $a₀       | ", abs(a - aᶠ) < abs(a - a₀) ? "    Yes!" : "   No!")
    println("     b     |   $b\t|    $bᶠ  \t|    $b₀       | ", abs(b - bᶠ) < abs(b - b₀) ? "    Yes!" : "   No!")
    println("     c     |   $c\t|    $cᶠ  \t|    $c₀       | ", abs(c - cᶠ) < abs(c - c₀) ? "    Yes!" : "   No!")
    println("     ρ°    |   $ρ\t|    $ρᶠ  \t|    $ρ₀       | ", abs(ρ - ρᶠ) < abs(ρ - ρ₀) ? "    Yes!" : "   No!")
    println("     λ°    |   $λ\t|    $λᶠ  \t|    $λ₀       | ", abs(λ - λᶠ) < abs(λ - λ₀) ? "    Yes!" : "   No!")
    println("     ν°    |   $ν\t|    $νᶠ  \t|    $ν₀       | ", abs(ν - νᶠ) < abs(ν - ν₀) ? "    Yes!" : "   No!")
    println("     βx    |   $βx\t|    $βxᶠ \t|    $βx₀       | ", abs(βx - βxᶠ) < abs(βx - βx₀) ? "    Yes!" : "   No!")
    println("     βy    |   $βy\t|    $βyᶠ \t|    $βy₀       | ", abs(βy - βyᶠ) < abs(βy - βy₀) ? "    Yes!" : "   No!")
    println("     βz    |   $βz\t|    $βzᶠ \t|    $βz₀       | ", abs(βz - βzᶠ) < abs(βz - βz₀) ? "    Yes!" : "   No!")
    println("__________________________________________________________________________")

    s_off, ζ_off, β_off = 0.2, 6, 2.0  # 2σ
    ##### PLOT ##### 
    ap = plot(s_est[1, :], title = "Scale Factor a", label = false); ap = hline!([a₀], ls = :dot, label = false); ap = hline!([a], ls = :dash, label = false, ylim = [a - s_off, a + s_off]);
    bp = plot(s_est[2, :], title = "Scale Factor b", label = false); bp = hline!([b₀], ls = :dot, label = false); bp = hline!([b], ls = :dash, label = false, ylim = [b - s_off, b + s_off]);
    cp = plot(s_est[3, :], title = "Scale Factor c", label = false); cp = hline!([c₀], ls = :dot, label = false); cp = hline!([c], ls = :dash, label = false, ylim = [c - s_off, c + s_off]);
    
    ρp = plot(ζ_est[1, :], title = "Non-Ortho Angles ρ", label = false); ρp = hline!([ρ₀], ls = :dot, label = false); ρp = hline!([ρ], ls = :dash, label = false, ylim = [ρ - ζ_off, ρ + ζ_off]);
    λp = plot(ζ_est[2, :], title = "Non-Ortho Angles λ", label = false); λp = hline!([λ₀], ls = :dot, label = false); λp = hline!([λ], ls = :dash, label = false, ylim = [λ - ζ_off, λ + ζ_off]);
    νp = plot(ζ_est[3, :], title = "Non-Ortho Angles ν", label = false); νp = hline!([ν₀], ls = :dot, label = false); νp = hline!([ν], ls = :dash, label = false, ylim = [ν - ζ_off, ν + ζ_off]);
    
    xp = plot(β_est[1, :], title = "Bias x", label = false); xp = hline!([βx₀], ls = :dot, label = false); xp = hline!([βx], ls = :dash, label = false, ylim = [βx - β_off, βx + β_off]);
    yp = plot(β_est[2, :], title = "Bias y", label = false); yp = hline!([βy₀], ls = :dot, label = false); yp = hline!([βy], ls = :dash, label = false, ylim = [βy - β_off, βy + β_off]);
    zp = plot(β_est[3, :], title = "Bias z", label = false); zp = hline!([βz₀], ls = :dot, label = false); zp = hline!([βz], ls = :dash, label = false, ylim = [βz - β_off, βz + β_off]);
    
    display(plot(ap, bp, cp))
    display(plot(ρp, λp, νp))
    display(plot(xp, yp, zp))
    
    return nothing
end;
magnetometer_calibration_report(results::NamedTuple) = magnetometer_calibration_report(results[:sat_truth], results[:sat_ests])


function diode_calibration_report(sat_true::SATELLITE, est_hist::Vector{SATELLITE{6, T}}) where {T} 
    N = size(est_hist, 1)
    C_est, α_est, ϵ_est = zeros(6, N), zeros(6, N), zeros(6, N)

    for i = 1:N
        C_est[:, i] = est_hist[i].diodes.calib_values
        α_est[:, i] = rad2deg.(est_hist[i].diodes.azi_angles)
        ϵ_est[:, i] = rad2deg.(est_hist[i].diodes.elev_angles)
    end
    
    C₀, α₀, ϵ₀ = est_hist[1].diodes.calib_values, rad2deg.(est_hist[1].diodes.azi_angles), rad2deg.(est_hist[1].diodes.elev_angles)
    Cf, αf, ϵf = est_hist[N].diodes.calib_values, rad2deg.(est_hist[N].diodes.azi_angles), rad2deg.(est_hist[N].diodes.elev_angles)
    C, α, ϵ = sat_true.diodes.calib_values, rad2deg.(sat_true.diodes.azi_angles), rad2deg.(sat_true.diodes.elev_angles)
    
    println("_____________________________________DIODES______________________________________")
    println("___DIODE___|_______Truth (C,α,ϵ)________|__________Final Guess__________|_________Init Guess__________|_____Improved?___")
    for i = 1:6 
        print("     $i     |    $(round(C[i], digits = 2)), $(round(α[i], digits = 2)), $(round(ϵ[i], digits = 2))  \t|")
        print("   $(round(Cf[i], digits = 2)), $(round(αf[i], digits = 2)), $(round(ϵf[i], digits = 2))    \t|")
        print("   $(round(C₀[i], digits = 2)), $(round(α₀[i], digits = 2)), $(round(ϵ₀[i], digits = 2))      \t|")
        print( (abs(C₀[i] - C[i]) ≤ abs(Cf[i] - C[i])) ? " No!, " : " Yes!, ")
        print( (abs(α₀[i] - α[i]) ≤ abs(αf[i] - α[i])) ? " No!, " : " Yes!, ")
        println( (abs(ϵ₀[i] - ϵ[i]) ≤ abs(ϵf[i] - ϵ[i])) ? " No! " : " Yes! ")
    end
    println("__________________________________________________________________________________")
    

    ##### PLOT #####
    C_off, ϵ_off, α_off = 0.25, 6.0, 6.0
    Cps, ϵps, αps = [], [], []
    for i = 1:6
        plot(C_est[i, :], title = "Scale Factor (C)", label = false); 
            hline!([C₀[i]], ls = :dot, label = false); Cp = hline!([C[i]], ls = :dash, label = false, ylim = [C[i] - C_off, C[i] + C_off])
        
        plot(ϵ_est[i, :], title = "Elevation Angle (ϵ)", label = false); 
            hline!([ϵ₀[i]], ls = :dot, label = false); ϵp = hline!([ϵ[i]], ls = :dash, label = false, ylim = [ϵ[i] - ϵ_off, ϵ[i] + ϵ_off])
        
        plot(α_est[i, :], title = "Azimuth Angle (α)", label = false); 
            hline!([α₀[i]], ls = :dot, label = false); αp = hline!([α[i]], ls = :dash, label = false, ylim = [α[i] - α_off, α[i] + α_off])
    

        push!(Cps, Cp)
        push!(ϵps, ϵp)
        push!(αps, αp)
    end

    # Adjust the labels 
    Cps[2].series_list[1][:label] = "Estimate"; Cps[2].series_list[2][:label] = "Initial Guess"; Cps[2].series_list[3][:label] = "Truth"
    ϵps[2].series_list[1][:label] = "Estimate"; ϵps[2].series_list[2][:label] = "Initial Guess"; ϵps[2].series_list[3][:label] = "Truth"
    αps[2].series_list[1][:label] = "Estimate"; αps[2].series_list[2][:label] = "Initial Guess"; αps[2].series_list[3][:label] = "Truth"

    display(plot(Cps..., layout = (3, 2)))
    display(plot(ϵps..., layout = (3, 2)))
    display(plot(αps..., layout = (3, 2)))
end;
diode_calibration_report(results::NamedTuple) = diode_calibration_report(results[:sat_truth], results[:sat_ests])

function mekf_report(states::Vector{STATE{T}}, est_hist::Vector{SATELLITE{6, T}}; verbose = true) where {T}
    N = size(states, 1)
    qs = [states[i].q for i = 1:N]; qs = reduce(hcat, qs)'; 
    βs = [states[i].β for i = 1:N]; βs = reduce(hcat, βs)';
    q̂s = [est_hist[i].state.q for i = 1:N]; q̂s = reduce(hcat, q̂s)';
    β̂s = [est_hist[i].state.β for i = 1:N]; β̂s = reduce(hcat, β̂s)';

    qErrs = [norm(cayley_map(qs[i, :], q̂s[i, :])) for i = 1:N];
    βErrs = [norm(βs[i, :] - β̂s[i, :]) for i = 1:N];

    if verbose
        plot( qs, title = "MEKF Report: q", c = [:red :orange :blue :green]);
        plot!(q̂s, c = [:red :orange :blue :green], ls = :dash, label = false);
        display( plot!(qErrs, label = false, c = :black, ls = :dot, ylim = (-1.5, 1.5), lw = 2) )

        plot( βs, title = "MEKF Report: β", c = [:red :blue :green]);
        plot!(β̂s, c = [:red :blue :green], ls = :dash, label = false);
        display( plot!(βErrs, label = false, c = :black, ls = :dot) )
    end
    # return qs, q̂s, βs, β̂s

    # Return average quaternion error of last ten guesses (?)
    return sum(qErrs[ (end-9):end ]) / 10
end;
mekf_report(results::NamedTuple, start = 1, stop = nothing; kwargs...) = (isnothing(stop)) ? mekf_report(results[:states][start:end] , results[:sat_ests][start:end]; kwargs...) :
                                                                                             mekf_report(results[:states][start:stop], results[:sat_ests][start:stop])

function detumbler_report(states, sensors; τ₁ = deg2rad(15), τ₂ = deg2rad(8) )
    N = size(states, 1)
    ωs = [states[i].ω for i = 1:N]; ωs = reduce(hcat, ωs)'
    ω̂s = [sensors[i].gyro for i = 1:N]; ω̂s = reduce(hcat, ω̂s)'
    nω = [norm(ωs[i, :]) for i = 1:N];
    nω̂ = [norm(ω̂s[i, :]) for i = 1:N];

    plot(ωs, c = [:red :blue :green], label = ["ωx" "ωy" "ωz"]); hline!([τ₁], c = :gray, ls = :dot, label = "τ₀"); 
        hline!([τ₂], c = :gray, ls = :dot, label = "τf"); 
        display(plot!(nω, ls = :dash, c = :black, label = "Mag", xlabel = "Step", ylabel = "Ang Velocity (rad/s)", title = "True Angular Velocity"))

    plot(ω̂s, c = [:red :blue :green], label = ["ωx" "ωy" "ωz"]); hline!([τ₁], c = :gray, ls = :dot, label = "τ₀"); 
        hline!([τ₂], c = :gray, ls = :dot, label = "τf"); 
        display(plot!(nω̂, ls = :dash, c = :black, label = "Mag", xlabel = "Step", ylabel = "Ang Velocity (rad/s)", title = "Measured Angular Velocity"))

    return ωs, ω̂s
end;
detumbler_report(results; kwargs...) = detumbler_report(results[:states], results[:sensors]; kwargs...)

function evaluate_diode_cal(sensors::Vector{SENSORS{6, T}}, truths::Vector{GROUND_TRUTH{6, T}}, d0::DIODES, df::DIODES; verbose = true) where {T}
    # Deal with numerical errors 
    r_acos(x) = (x ≈  1) ? zero(x)    : 
                (x ≈ -1) ? one(x) * π : acos(x)

    no_eclipse = [norm(sensors[i].diodes) > 0.02 for i = 1:size(sensors, 1)]
    
    
    sensors = sensors[ no_eclipse ]; # Ignore eclipse
    truths  = truths[  no_eclipse ]; 

    N = size(truths, 1)

    ŝ0 = [estimate_sun_vector(sensors[i], d0) for i = 1:N];
    ŝf = [estimate_sun_vector(sensors[i], df) for i = 1:N];
    sᴮ = [truths[i].ŝᴮ for i = 1:N];

    e0 = [ rad2deg( r_acos(ŝ0[i]' * sᴮ[i])) for i = 1:N]
    ef = [ rad2deg( r_acos(ŝf[i]' * sᴮ[i])) for i = 1:N]
    # e0 = [norm(sᴮ[i] - ŝ0[i]) for i = 1:N]
    # ef = [norm(sᴮ[i] - ŝf[i]) for i = 1:N]

    # Remove NaNs from eclipses
    e0 = e0[.!isnan.(e0)]
    ef = ef[.!isnan.(ef)]

    μ0 = round(sum(e0) / size(e0, 1), digits = 3); σ0 = round(std(e0), digits = 3)
    μf = round(sum(ef) / size(ef, 1), digits = 3); σf = round(std(ef), digits = 3)

    if verbose
        histogram(e0, title = "Initial Error (sᴮ)", xlabel = "Err (deg)", ylabel = "Frequency", label = false);
        hline!([0 ], label = "σ0 = $σ0");
        vline!([μ0], label = "μ0 = $μ0");
        pe0 = vline!([μf], ls = :dash, label = "(μf = $μf)");

        histogram(ef, title = "Final Error (sᴮ)", xlabel = "Err (deg)", ylabel = "Frequency", label = false);
        hline!([0 ], label = "σf = $σf");
        vline!([μf], label = "μf = $μf");
        pef = vline!([μ0], ls = :dash, label = "(μ0 = $μ0)");

        display(plot(pe0, pef, layout = (2, 1)))

        println("----- DIODE CALIBRATION -----")
        println("Initial Error (deg): μ: $μ0,  σ: $σ0")
        println("Final   Error (deg): μ: $μf,  σ: $σf")
        println()
    end

    return μ0, μf
end
evaluate_diode_cal(results::NamedTuple; kwargs...) = evaluate_diode_cal(results[:sensors], results[:truths], results[:sat_ests][1].diodes, results[:sat_ests][end].diodes; kwargs...)


# CAN ONLY BE RUN 𝑏𝑒𝑓𝑜𝜖𝑒 it has been calibrated, else it is already corrected 
function evaluate_mag_cal(sensors::Vector{SENSORS{6, T}}, truths::Vector{GROUND_TRUTH{6, T}}, modes::Vector{Operation_mode}, 
                            sat0::SATELLITE, satf::SATELLITE; verbose = true, force_show = false) where {T}
    # Deal with numerical errors 
    r_acos(x) = (x ≈  1) ? zero(x)    : 
                (x ≈ -1) ? one(x) * π : acos(x)
    
    N = size(truths, 1)

    B̂0 = [correct_magnetometer(sat0, sensors[i].magnetometer) for i = 1:N]
    B̂f = [correct_magnetometer(satf, sensors[i].magnetometer) for i = 1:N]
    Bᴮ = [truths[i].Bᴮ for i = 1:N]

    e0 = [ rad2deg( r_acos(normalize(B̂0[i])' * normalize(Bᴮ[i]))) for i = 1:N]
    ef = [ rad2deg( r_acos(normalize(B̂f[i])' * normalize(Bᴮ[i]))) for i = 1:N]

    # Once the magnetometer has been calibrated, we have nothing to compare it against (?)
    if !(force_show)
        e0 = e0[modes .== mag_cal]
        ef = ef[modes .== mag_cal]
    end
    N  = size(e0, 1)

    if N > 2
        μ0 = sum(e0) / N; σ0 = std(e0)
        μf = sum(ef) / N; σf = std(ef)

        if verbose
            println("----- MAGNETOMETER CALIBRATION -----")
            println("Initial Error (deg): μ: $μ0,  σ: $σ0")
            println("Final   Error (deg): μ: $μf,  σ: $σf")
            println()
        end

        return μ0, μf
    else 
        println("Error - Accuracy of Mag vec in Body frame cannot be evaluated!")
        return 1000, 1000
    end
end
evaluate_mag_cal(results::NamedTuple; kwargs...) = evaluate_mag_cal(results[:sensors], results[:truths], results[:modes], results[:sat_ests][1], results[:sat_ests][end]; kwargs...)

# @warn "Monte carlo not using albedo (for speed)! And testing mekf only rn"

function monte_carlo(N = 10)
    # Average pointing error, sun error, mag error (?)
    eqs, ess, eBs = zeros(N), zeros(N), zeros(N)
    es₀s, eB₀s    = zeros(N), zeros(N)
    for i = 1:N
        println("\n---------  $i --------")
        # results = main(; verbose = false, num_orbits = 1.5, initial_state = diode_cal, σβ = 0.0, σB = 0.0, σ_gyro = 0.0, σr = 0.0, σ_current = 0.0);
        results = main(; initial_mode = mekf, num_orbits = 1.0, verbose = false, use_albedo = true)
        eq = mekf_report(results; verbose = false)
        es₀, es = evaluate_diode_cal(results; verbose = false)
        eB₀, eB = evaluate_mag_cal(results; verbose = false, force_show = true)
        eqs[i] = eq 
        ess[i] = es 
        eBs[i] = eB
        es₀s[i] = es₀
        eB₀s[i] = eB₀
    end

    sDiff = (ess - es₀s);
    bDiff = (eBs - eB₀s);
    μsf, σsf = round(mean(ess), digits = 2),   round(std(ess), digits = 2);
    μs0, σs0 = round(mean(es₀s), digits = 2),  round(std(es₀s), digits = 2);
    μsd, σsd = round(mean(sDiff), digits = 2), round(std(sDiff), digits = 2);
    μbf, σbf = round(mean(eBs), digits = 2),   round(std(eBs), digits = 2);
    μb0, σb0 = round(mean(eB₀s), digits = 2),  round(std(eB₀s), digits = 2);
    μbd, σbd = round(mean(bDiff), digits = 2), round(std(bDiff), digits = 2);
    μϕ, σϕ   = round(mean(eqs), digits = 3), round(std(eqs), digits = 3);

    header  = (["Vector", "Initial", "Final", "Difference (fin - init)"], ["N = $N", "μ° (σ°)", "μ° (σ°)", "μ (σ)"]);
    vectors = ["Sun", "Mag", "Attitude"];
    inits   = ["$μs0 ($σs0)", "$μb0 ($σb0)", "-----"];
    finals  = ["$μsf ($σsf)", "$μbf ($σbf)", "$μϕ ($σϕ)"];
    diffs   = ["$μsd ($σsd)", "$μbd ($σbd)", "-----"];
    data    = hcat(vectors, inits, finals, diffs);
    display(pretty_table(data; header = header, header_crayon = crayon"yellow bold"));

    return eqs, ess, eBs, es₀s, eB₀s
end



# NOTE that we are using a sqrt KF, so ``Covariance'' is already a stdev 
function mag_self_consistency(results, t::Symbol = :a; start = 1, stop = nothing, ds = 1, kwargs...)
    sat_true, est_hist = results[:sat_truth], results[:sat_ests];
    N = size(est_hist, 1)

    function wrap(diff; radians = true)

        offset = (radians) ? pi : 180
        while diff > offset 
            diff -= 2 * offset 
        end

        while diff < -offset 
            diff += 2 * offset 
        end

        return diff
    end

    s_err, ζ_err, β_err = zeros(N, 3), zeros(N, 3), zeros(N, 3);
    σs, σζ, σβ = zeros(N, 3), zeros(N, 3), zeros(N, 3);

    a, b, c    = sat_true.magnetometer.scale_factors;
    ρ, λ, ν    = sat_true.magnetometer.non_ortho_angles;
    βx, βy, βz = sat_true.magnetometer.bias;
    
    for i = 1:N
        s_err[i, :] = [a, b, c] .- est_hist[i].magnetometer.scale_factors
        ζ_err[i, :] = wrap.([ρ, λ, ν] .- est_hist[i].magnetometer.non_ortho_angles)
        β_err[i, :] = [βx, βy, βz] .- est_hist[i].magnetometer.bias

        Σchol = est_hist[i].covariance 
        Σ = Σchol' * Σchol 
        σs[i, :] .= sqrt.(diag(Σ[7:9, 7:9]))
        σζ[i, :] .= sqrt.(diag(Σ[10:12, 10:12]))
        σβ[i, :] .= sqrt.(diag(Σ[13:15, 13:15]))

        # σs[i, :] = ([est_hist[i].covariance[7,   7], est_hist[i].covariance[8,   8], est_hist[i].covariance[9,   9]])
        # σζ[i, :] = ([est_hist[i].covariance[10, 10], est_hist[i].covariance[11, 11], est_hist[i].covariance[12, 12]])
        # σβ[i, :] = ([est_hist[i].covariance[13, 13], est_hist[i].covariance[14, 14], est_hist[i].covariance[15, 15]])
    end

    plot(s_err, title = "Error - Scale Factors", label = ["a" "b" "c"], xlabel = "Index", ylabel = "Value", layout = 3);
    plot!(      3 * σs, c = :red, layout = 3, label = "3σ");
    ps = plot!(-3 * σs, c = :red, layout = 3, label = false)

    plot(rad2deg.(ζ_err), title = "Error - Non-Ortho Angles", label = ["ρ" "λ" "ν"], xlabel = "Index", ylabel = "Value", layout = 3);
    plot!(      3 * rad2deg.(σζ), c = :red, layout = 3, label = "3σ");
    pζ = plot!(-3 * rad2deg.(σζ), c = :red, layout = 3, label = false)

    plot(β_err, title = "Error - Bias", label = ["x" "y" "z"], xlabel = "Index", ylabel = "Value", layout = 3);
    plot!(      3 * σβ, c = :red, layout = 3, label = "3σ");
    pβ = plot!(-3 * σβ, c = :red, layout = 3, label = false)

    display(plot(ps))
    display(plot(pζ))
    display(plot(pβ))

    return nothing
end

function diode_self_consistency(results, t::Symbol = :a; start = 1, stop = nothing, ds = 1, kwargs...)
    sat_true, est_hist = results[:sat_truth], results[:sat_ests];
    N = size(est_hist, 1)
    Nd = 6

    function wrap(diff; radians = true)

        offset = (radians) ? pi : 180
        while diff > offset 
            diff -= 2 * offset 
        end

        while diff < -offset 
            diff += 2 * offset 
        end

        return diff
    end

    
    C_err, α_err, ϵ_err = zeros(N, Nd), zeros(N, Nd), zeros(N, Nd);
    σC, σα, σϵ = zeros(N, Nd), zeros(N, Nd), zeros(N, Nd);
    Cs, αs, ϵs = sat_true.diodes.calib_values, sat_true.diodes.azi_angles, sat_true.diodes.elev_angles;

    for i = 1:N
        C_err[i, :] .= Cs .- est_hist[i].diodes.calib_values;
        α_err[i, :] .= wrap.( αs - est_hist[i].diodes.azi_angles);
        ϵ_err[i, :] .= wrap.( ϵs - est_hist[i].diodes.elev_angles);

        # NO sqrt. because it is already a square root KF (but keep abs)
        Σchol = est_hist[i].covariance 
        Σ = Σchol' * Σchol 

        σC[i, :] .= sqrt.(diag(Σ[16:21, 16:21]))  # abs.(diag(est_hist[i].covariance[16:21, 16:21]));
        σα[i, :] .= sqrt.(diag(Σ[22:27, 22:27]))  #abs.(diag(est_hist[i].covariance[22:27, 22:27]));
        σϵ[i, :] .= sqrt.(diag(Σ[28:33, 28:33]))  #abs.(diag(est_hist[i].covariance[28:33, 28:33]));
    end

    plot(C_err, title = "Error - Scale Factors", xlabel = "Index", ylabel = "Value", layout = 6);
    plot!(      3 * σC, c = :red, layout = 6, label = "3σ");
    pC = plot!(-3 * σC, c = :red, layout = 6, label = false)

    plot(rad2deg.(α_err), title = "Error - Azimuth Angles",   xlabel = "Index", ylabel = "Value", layout = 6);
    plot!(      3 * rad2deg.(σα), c = :red, layout = 6, label = "3σ");
    pα = plot!(-3 * rad2deg.(σα), c = :red, layout = 6, label = false)

    plot(rad2deg.(ϵ_err), title = "Error - Elevation Angles", xlabel = "Index", ylabel = "Value", layout = 6);
    plot!(      3 * rad2deg.(σϵ), c = :red, layout = 6, label = "3σ");
    pϵ = plot!(-3 * rad2deg.(σϵ), c = :red, layout = 6, label = false)

    display(plot(pC))
    display(plot(pα))
    display(plot(pϵ))

    return nothing
end

function mekf_self_consistency(results)
    # Deal with numerical errors 
    r_acos(x) = (x ≈  1) ? zero(x)    : 
                (x ≈ -1) ? one(x) * π : acos(x)

                
    x = results[:states]
    N = size(x, 1)
    x̂ = results[:sat_ests] # [results[:sat_ests].state for i = 1:N]

    ϕ_err, β_err = zeros(N), zeros(N, 3);
    σϕ, σβ = zeros(N), zeros(N, 3);
    vᴮ = normalize([1, 1, 1]);
    for i = 1:N
        Q, Q̂ = quat2rot(x[i].q), quat2rot(x̂[i].state.q)
        ϕ_err[i] = rad2deg(r_acos((Q * vᴮ)' * (Q̂ * vᴮ)))
        β_err[i, :] .= x[i].β - x̂[i].state.β

        Σchol = x̂[i].covariance
        Σ = Σchol' * Σchol
         
        σϕ[i] = rad2deg(norm(sqrt.(diag(Σ[1:3, 1:3]))))  #  rad2deg(norm(diag(x̂[i].covariance[1:3, 1:3])))
        σβ[i, :] .= sqrt.(diag(Σ[4:6, 4:6]))  #  abs.(diag(x̂[i].covariance[4:6, 4:6]))
    end

    plot(ϕ_err, title = "Error - Attitude");
    plot!(   3 * σϕ);
    pϕ = plot!(-3 * σϕ)

    plot(β_err, title = "Error - Gyro Bias", layout = 3);
    plot!(      3 * σβ, c = :red, layout = 3, label = "3 * σ");
    pβ = plot!(-3 * σβ, c = :red, layout = 3, label = false)



end



# # Get stdev of q est 
# @testset "σϕ" begin
#     alb = get_albedo(2)
#     dt  = 0.2
#     vᴮ  = [1, 1, 1]; vᴮ /= norm(vᴮ)
#     errs = []

#     for i = 1:15000
#         t    = Epoch(2020, 1, 1) + (60 * 60 * 24) * rand(1:365)
#         sat  = SATELLITE()
#         sat_est = SATELLITE(; ideal = true)
#         x, _ = get_initial_state()

#         truth, sensors, ecl, noise = generate_measurements(sat, alb, x, t, dt; use_albedo = true);
#         if ecl > 0.1
#             q = run_triad(sensors, sat_est, t, true; sᴮ_true = truth.ŝᴮ)  # x.q 
#             Q̂ = quat2rot(q); Q = quat2rot(x.q)
#             v̂ᴵ = Q̂ * vᴮ; vᴵ = Q * vᴮ;
#             er = rad2deg(acos(v̂ᴵ' * vᴵ))
#             push!(errs, er)
#         end
#     end
# end


function figs_for_poster(results)
    r_acos(x) = (x ≈  1) ? zero(x)    : 
                (x ≈ -1) ? one(x) * π : acos(x)

    N = size(results[:states], 1);
    sensors = results[:sensors]; 
    truths  = results[:truths]; 
    ests    = results[:sat_ests];
    states  = results[:states];

    #### Diode Cal 
    ŝs = [estimate_sun_vector(sensors[i], ests[i].diodes) for i = 1:N];
    sᴮ = [truths[i].ŝᴮ for i = 1:N];
    es = [ rad2deg( r_acos(ŝs[i]' * sᴮ[i])) for i = 1:5:N];

    idx = 1
    Cs = [results[:sat_ests][i].diodes.calib_values[idx] for i = 1:5:N];         # Down sample to make it in seconds
    αs = [rad2deg(results[:sat_ests][i].diodes.azi_angles[idx])  for i = 1:5:N];
    ϵs = [rad2deg(results[:sat_ests][i].diodes.elev_angles[idx]) for i = 1:5:N];

    C₀, α₀, ϵ₀ = Cs[idx], αs[idx], ϵs[idx]
    C, α, ϵ = results[:sat_truth].diodes.calib_values[idx], rad2deg(results[:sat_truth].diodes.azi_angles[idx]), rad2deg(results[:sat_truth].diodes.elev_angles[idx]);

    # pc = plot(Cs, title = "Calibration Value", label = "Estimate", ylabel = "Magnitude", xlabel = "Time (s)");  hline!([C₀], ls = :dot, c = :red, label = "Initial Guess"); hline!([C], ls = :dash, c = :green, label = "True Value", legend = false);
    # pa = plot(αs, title = "Azimuth Angle", label = "Estimate", ylabel = "Angle (deg)", xlabel = "Time (s)");    hline!([α₀], ls = :dot, c = :red, label = "Initial Guess"); hline!([α], ls = :dash, c = :green, label = "True Value");
    # pe = plot(ϵs, title = "Elevation Angle", label = "Estimate", ylabel = "Angle (deg)", xlabel = "Time (s)");  hline!([ϵ₀], ls = :dot, c = :red, label = "Initial Guess"); hline!([ϵ], ls = :dash, c = :green, label = "True Value", legend = false);
    # ps = plot(es, title = "Estimation Error", label = "Sun Vector", ylabel = "Error (deg)", xlabel = "Time (s)")

    pc = plot(Cs, label = false); hline!([C₀], ls = :dot, c = :red, label = false); hline!([C], ls = :dash, c = :green, label = false, legend = false);
    pa = plot(αs, label = false); hline!([α₀], ls = :dot, c = :red, label = false); hline!([α], ls = :dash, c = :green, label = false);
    pe = plot(ϵs, label = false); hline!([ϵ₀], ls = :dot, c = :red, label = false); hline!([ϵ], ls = :dash, c = :green, label = false, legend = false);
    # ps = plot(es, label = false); 
    cal_plot = plot(pc, pa, pe, xlim = [1, 2500], legend = false)

    ### Magnetometer Calibration 
    B̂0 = [correct_magnetometer(ests[1],   sensors[i].magnetometer) for i = 1:N];
    B̂f = [correct_magnetometer(ests[end], sensors[i].magnetometer) for i = 1:N];
    Bᴮ = [truths[i].Bᴮ for i = 1:N];

    e0 = [ rad2deg( r_acos(normalize(B̂0[i])' * normalize(Bᴮ[i]))) for i = 1:N];
    ef = [ rad2deg( r_acos(normalize(B̂f[i])' * normalize(Bᴮ[i]))) for i = 1:N];

    μ0 = round(sum(e0) / N, digits = 2); μf = round(sum(ef) / N, digits = 2);
    h0 = histogram(e0, normalize = true, ylabel = "Frequency", xlabel = "Error", title = "Before Calibration", label = "μ = $μ0");
    hf = histogram(ef, normalize = true, ylabel = "Frequency", xlabel = "Error", title = "After Calibration", label = "μ = $μf");
    mag_plot = plot(h0, hf, plot_title = "Magnetometer Calibration", layout = (2, 1))


    ### MEKF 
    qs = [states[i].q for i = 1:5:N]; qs = reduce(hcat, qs)'; 
    q̂s = [est_hist[i].state.q for i = 1:5:N]; q̂s = reduce(hcat, q̂s)';
    mekf_plot = plot(qErrs, label = false, title = "Attitude Error", xlabel = "Time (s)", ylabel = "Error magnitude", xlim = [1, 2500])


    ### Detumbling 
    modes = results[:modes]
    Nd = findall( (modes .== detumble) .== 0)[1]
    ωs = [states[i].ω for i = 1:Nd]; ωs = reduce(hcat, ωs)';
    ω̂s = [sensors[i].gyro for i = 1:Nd]; ω̂s = reduce(hcat, ω̂s)';
    nω = [norm(ωs[i, :]) for i = 1:5:Nd];
    nω̂ = [norm(ω̂s[i, :]) for i = 1:5:Nd];
    τ₁ = 15 # deg2rad(15)
    τ₂ = 5
    
    plot(rad2deg.(nω̂ ), label = "Measured", title = "Gyroscope");
    plot!(rad2deg.(nω), label = "Truth", xlabel = "Time (s)", ylabel = "Magnitude (deg/s)", ls = :dash);
    hline!([τ₁], label = "Initial Detumble", ls = :dot, c = :gray)
    det_plot = hline!([τ₂], label = "Final Detumble", ls = :dot, c = :gray)


    return cal_plot, mag_plot, mekf_plot, det_plot
end