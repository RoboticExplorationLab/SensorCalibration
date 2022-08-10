# [src/reports.jl]

"""
    Contains functions for evaluating the performance of the Sim/Filter setup 
"""

function magnetometer_calibration_report(sat_true, sat_est, sat_init) 

    aᶠ, bᶠ, cᶠ = round.(sat_est.magnetometer.scale_factors, sigdigits = 3)
    ρᶠ, λᶠ, ϕᶠ = round.(sat_est.magnetometer.non_ortho_angles, sigdigits = 3)
    βxᶠ, βyᶠ, βzᶠ = round.(sat_est.magnetometer.bias, sigdigits = 3)

    a, b, c = round.(sat_true.magnetometer.scale_factors, sigdigits = 3)
    ρ, λ, ϕ = round.(sat_true.magnetometer.non_ortho_angles, sigdigits = 3)
    βx, βy, βz = round.(sat_true.magnetometer.bias, sigdigits = 3)

    a₀, b₀, c₀ = round.(sat_init.magnetometer.scale_factors, sigdigits = 3)
    ρ₀, λ₀, ϕ₀ = round.(sat_init.magnetometer.non_ortho_angles, sigdigits = 3)
    βx₀, βy₀, βz₀ = round.(sat_init.magnetometer.bias, sigdigits = 3)

    println("__________________________________________________________________________")
    println("___PARAM___|___Truth____|__Final Guess__|__Init Guess__|__Improved?__")
    println("     a     |   $a\t|    $aᶠ  \t|    $a₀       | ", abs(a - aᶠ) < abs(a - a₀) ? "    Yes!" : "   No!")
    println("     b     |   $b\t|    $bᶠ  \t|    $b₀       | ", abs(b - bᶠ) < abs(b - b₀) ? "    Yes!" : "   No!")
    println("     c     |   $c\t|    $cᶠ  \t|    $c₀       | ", abs(c - cᶠ) < abs(c - c₀) ? "    Yes!" : "   No!")
    println("     ρ°    |   $ρ\t|    $ρᶠ  \t|    $ρ₀       | ", abs(ρ - ρᶠ) < abs(ρ - ρ₀) ? "    Yes!" : "   No!")
    println("     λ°    |   $λ\t|    $λᶠ  \t|    $λ₀       | ", abs(λ - λᶠ) < abs(λ - λ₀) ? "    Yes!" : "   No!")
    println("     ϕ°    |   $ϕ\t|    $ϕᶠ  \t|    $ϕ₀       | ", abs(ϕ - ϕᶠ) < abs(ϕ - ϕ₀) ? "    Yes!" : "   No!")
    println("     βx    |   $βx\t|    $βxᶠ \t|    $βx₀       | ", abs(βx - βxᶠ) < abs(βx - βx₀) ? "    Yes!" : "   No!")
    println("     βy    |   $βy\t|    $βyᶠ \t|    $βy₀       | ", abs(βy - βyᶠ) < abs(βy - βy₀) ? "    Yes!" : "   No!")
    println("     βz    |   $βz\t|    $βzᶠ \t|    $βz₀       | ", abs(βz - βzᶠ) < abs(βz - βz₀) ? "    Yes!" : "   No!")
    println("__________________________________________________________________________")

    return nothing
end;
magnetometer_calibration_report(results::NamedTuple) = magnetometer_calibration_report(results[:sat_truth], results[:sat_ests][end], results[:sat_ests][1])


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
    C_off, ϵ_off, α_off = 0.2, 6.0, 6.0
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
    qs = [states[i].q for i = 1:5:N]; qs = reduce(hcat, qs)'; 
    βs = [states[i].β for i = 1:5:N]; βs = reduce(hcat, βs)';
    q̂s = [est_hist[i].state.q for i = 1:5:N]; q̂s = reduce(hcat, q̂s)';
    β̂s = [est_hist[i].state.β for i = 1:5:N]; β̂s = reduce(hcat, β̂s)';
    N = size(qs, 1)

    # Cayley map ≈ 1/2 of the error, in radians, as long as it is small (~10-15 deg)
    qErrs = [rad2deg(2 * norm(cayley_map(qs[i, :], q̂s[i, :]))) for i = 1:N];
    βErrs = [norm(βs[i, :] - β̂s[i, :]) for i = 1:N];

    if maximum(qErrs) > 15
        @warn "Quaternion error is large, Cayley map may not be very accurate"
    end

    if verbose
        plot( qs, title = "MEKF Report: q", c = [:red :orange :blue :green]);
        plot!(q̂s, c = [:red :orange :blue :green], ls = :dash, label = false);
        display( plot!(qErrs, label = false, ylabel = "Error (deg)", c = :black, ls = :dot, ylim = (-1.5, 1.5), lw = 2) )

        plot( βs, title = "MEKF Report: β", c = [:red :blue :green]);
        plot!(β̂s, c = [:red :blue :green], ls = :dash, label = false);
        display( plot!(βErrs, xlabel = "Time (s)", ylabel = "Bias (rad/s)", label = "Error", c = :black, ls = :dot) )
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

    no_eclipse = [norm(sensors[i].diodes) > 0.1 for i = 1:size(sensors, 1)]

    sensors = sensors[ no_eclipse ]; # Ignore eclipse
    truths  = truths[  no_eclipse ]; 
    
    N = size(truths, 1)

    ŝ0 = [estimate_sun_vector(sensors[i], d0) for i = 1:N]
    ŝf = [estimate_sun_vector(sensors[i], df) for i = 1:N]
    sᴮ = [truths[i].ŝᴮ for i = 1:N]

    e0 = [ rad2deg( acos(ŝ0[i]' * sᴮ[i])) for i = 1:N]
    ef = [ rad2deg( acos(ŝf[i]' * sᴮ[i])) for i = 1:N]
    # e0 = [norm(sᴮ[i] - ŝ0[i]) for i = 1:N]
    # ef = [norm(sᴮ[i] - ŝf[i]) for i = 1:N]

    # Remove NaNs from eclipses
    e0 = e0[.!isnan.(e0)]
    ef = ef[.!isnan.(ef)]

    μ0 = round(sum(e0) / size(e0, 1), digits = 3); σ0 = round(std(e0), digits = 3)
    μf = round(sum(ef) / size(ef, 1), digits = 3); σf = round(std(ef), digits = 3)

    if verbose
        pe0 = histogram(e0, title = "Initial Error (sᴮ)", bins = 0:20, normalize = true, xlabel = "Err (deg)", ylabel = "Frequency", label = "μ = $(μ0)°");
        hline!([0 ], label = "σ0 = $σ0");
        vline!([μ0], label = "μ0 = $μ0");
        pe0 = vline!([μf], ls = :dash, label = "(μf = $μf)");

        pef = histogram(ef, title = "Final Error (sᴮ)", bins = 0:20, normalize = true, xlabel = "Err (deg)", ylabel = "Frequency", label = "μ = $(μf)°");
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
                            sat0::SATELLITE, satf::SATELLITE; verbose = true) where {T}
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
    e0 = e0[modes .== mag_cal]
    ef = ef[modes .== mag_cal]
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


function monte_carlo(N = 10)
    # Average pointing error, sun error, mag error (?)
    eqs, ess, eBs = zeros(N), zeros(N), zeros(N)
    es₀s, eB₀s    = zeros(N), zeros(N)
    for i = 1:N
        println("\n---------  $i --------")
        results = main(; initial_mode = mekf, verbose = false, num_orbits = 1.0)  # mekf + 1 orbit + force mag_cal

        eq = mekf_report(results; verbose = false)
        es₀, es = evaluate_diode_cal(results; verbose = false)
        eB₀, eB = evaluate_mag_cal(results; verbose = false)
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

    header  = (["Vector", "Initial", "Final", "Difference (fin - init)"], ["N = $N", "μ° (σ°)", "μ° (σ°)", "μ° (σ°)"]);
    vectors = ["Sun", "Mag", "Attitude"];
    inits   = ["$μs0 ($σs0)", "$μb0 ($σb0)", "-----"];
    finals  = ["$μsf ($σsf)", "$μbf ($σbf)", "$μϕ ($σϕ)"];
    diffs   = ["$μsd ($σsd)", "$μbd ($σbd)", "-----"];
    data    = hcat(vectors, inits, finals, diffs);

    @show data
    display(pretty_table(data))
    display(pretty_table(data; header = header, header_crayon = crayon"yellow bold"));
    

    @infiltrate
    return eqs, ess, eBs, es₀s, eB₀s
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

        σC[i, :] .= sqrt.(diag(Σ[7:12, 7:12]))  #abs.(diag(est_hist[i].covariance[16:21, 16:21]));
        σα[i, :] .= sqrt.(diag(Σ[13:18, 13:18]))  #abs.(diag(est_hist[i].covariance[22:27, 22:27]));
        σϵ[i, :] .= sqrt.(diag(Σ[19:24, 19:24]))  #abs.(diag(est_hist[i].covariance[28:33, 28:33]));
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

function cal_gif(sat_true::SATELLITE, est_hist::Vector{SATELLITE{6, T}}) where {T} 
    N = Int(floor(size(est_hist, 1) / 5))
    C_est, α_est, ϵ_est = zeros(6, N), zeros(6, N), zeros(6, N)

    for i = 1:N
        C_est[:, i] = est_hist[5 * i + 4].diodes.calib_values
        α_est[:, i] = rad2deg.(est_hist[5 * i + 4].diodes.azi_angles)
        ϵ_est[:, i] = rad2deg.(est_hist[5 * i + 4].diodes.elev_angles)
    end
    
    C₀, α₀, ϵ₀ = est_hist[1].diodes.calib_values, rad2deg.(est_hist[1].diodes.azi_angles), rad2deg.(est_hist[1].diodes.elev_angles)
    Cf, αf, ϵf = est_hist[N].diodes.calib_values, rad2deg.(est_hist[N].diodes.azi_angles), rad2deg.(est_hist[N].diodes.elev_angles)
    C, α, ϵ = sat_true.diodes.calib_values, rad2deg.(sat_true.diodes.azi_angles), rad2deg.(sat_true.diodes.elev_angles)

    idx = 1
    p1 = hline([C₀[idx]], title = "Scale Factor",  c = :red, label = false, ls = :dot); p1 = hline!([C[idx]], c = :green, label = false, ls = :dash);
    p1 = @gif for i = 1:400 #N
         plot!(C_est[idx, 1:i], xlim = [0, 400], label = false, c = :blue, ylim = [C[idx] - C_off, C[idx] + C_off])
        # p2 = plot!(α_est[idx, 1:i], xlim = [0, 400], label = false, c = :blue, ylim = [α[idx] - α_off, α[idx] + α_off])
    end every 10;

    p2 = hline([α₀[idx]], title = "Azimuth Angle", c = :red, label = false, ls = :dot); p2 = hline!([α[idx]], c = :green, label = false, ls = :dash);
    p2 = @gif for i = 1:400 #N
        # plot!(C_est[idx, 1:i], xlim = [0, 400], label = false, c = :blue, ylim = [C[idx] - C_off, C[idx] + C_off])
       plot!(α_est[idx, 1:i], xlim = [0, 400], label = false, c = :blue, ylim = [α[idx] - α_off, α[idx] + α_off])
    end every 10;

    p3 = hline([α₀[idx]], title = "Elevation Angle", c = :red, label = false, ls = :dot); p3 = hline!([ϵ[idx]], c = :green, label = false, ls = :dash);
    p3 = @gif for i = 1:400 #N
        # plot!(C_est[idx, 1:i], xlim = [0, 400], label = false, c = :blue, ylim = [C[idx] - C_off, C[idx] + C_off])
       plot!(ϵ_est[idx, 1:i], xlim = [0, 400], label = false, c = :blue, ylim = [ϵ[idx] - ϵ_off, ϵ[idx] + ϵ_off])
    end every 10;

    return p1, p2, p3
end;