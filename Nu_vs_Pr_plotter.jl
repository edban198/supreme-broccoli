using DataFrames, CairoMakie, CSV, GLM, StatsModels, Printf, LaTeXStrings, Colors

# --- Load data ---
df = CSV.read("data/Nu_vs_Pr.csv", DataFrame)
df = dropmissing(df)

# --- Get unique χ values ---
χ_values = sort(unique(df.chi))

# --- Plot setup ---
fig = Figure(size=(1000, 600))
ax = Axis(fig[1,1],
    xlabel = "Pr",
    ylabel = "Nu",
    xscale = log10,
    yscale = log10,
    title = "Nu vs Pr (log-log) for varying χ"
)

# --- Color palette ---
palette = [
    :dodgerblue,
    :crimson,
    :forestgreen,
    :darkorange,
    :mediumvioletred,
    :teal,
    :indigo
][1:length(χ_values)]  # Trim if fewer χ


# --- Loop over each χ ---
for (i, χ) in enumerate(χ_values)
    dfχ = sort(filter(:chi => ==(χ), df), :Prandtl)
    color = palette[i]

    # Scatter points
    scatter!(ax, dfχ.Prandtl, dfχ.Nu, label="χ = $χ", color=color)

    # Log-log fit
    dfχ_log = DataFrame(logPr = log10.(dfχ.Prandtl), logNu = log10.(dfχ.Nu))
    fit = lm(@formula(logNu ~ logPr), dfχ_log)
    slope, intercept = coef(fit)[2], coef(fit)[1]
    A = 10^intercept

    # Store equation as LaTeXString
    eqn_str = @sprintf("\\chi = %.0f:\\quad Nu \\approx %.3f \\cdot Pr^{%.3f}", χ, A, slope)
    eqn = LaTeXString(eqn_str)

    # Prediction line
    Pr_range = range(extrema(dfχ.Prandtl)..., length=200)
    Nu_fit = A .* Pr_range .^ slope
    lines!(ax, Pr_range, Nu_fit, color=color, linestyle=:dash)

    # Annotate on plot (space them out vertically)
    y_pos = 10 ^ (log10(maximum(dfχ.Nu)) - 0.2 * i)
    text!(ax, minimum(dfχ.Prandtl) * 1.1, y_pos, text=eqn, color=color, align=(:left, :center))

    @info "χ = $χ → Nu ≈ $(round(A, sigdigits=3)) * Pr^$(round(slope, digits=3))"
end

# Legend
axislegend(ax, position=:rb)

# Save
save("OUTPUTS/Nu_vs_Pr.png", fig)