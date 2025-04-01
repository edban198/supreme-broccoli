using CairoMakie
using CSV, DataFrames
using Statistics, GLM

# Data
df = CSV.read("data/flux_wind_forcing.csv", DataFrame)

F = df.flux_wind_forcing
Nu = df.Nu
γ = df.gamma
slurm_id = df.SLURM_id

# Separate by R/R_c
F1 = [F[i] for i in eachindex(γ) if γ[i] == 2.0]
Nu1 = [Nu[i] for i in eachindex(γ) if γ[i] == 2.0]

F2 = [F[i] for i in eachindex(γ) if γ[i] == 2.5]
Nu2 = [Nu[i] for i in eachindex(γ) if γ[i] == 2.5]

# Plot
fig = Figure(size=(800, 800))
ax = Axis(fig[1, 1],
          xlabel="Wind forcing (flux)",
          ylabel="Nu",
          title="Nu vs flux Wind Forcing at Different R/R_c"
)

scatter!(ax, F1, Nu1, label="R/R_c = 2")
scatter!(ax, F2, Nu2, label="R/R_c = 2.5")

axislegend(ax)

# Save figure
save("nu_vs_flux_wind_forcing.png", fig)

df_c = CSV.read("data/constant_wind_forcing.csv", DataFrame)

F_c = df_c.const_wind_forcing
Nu_c = df_c.Nu
γ_c = df_c.gamma
slurm_id_c = df_c.SLURM_id

F1_c = [F_c[i] for i in eachindex(γ_c) if γ_c[i] == 2.0]
Nu1_c = [Nu_c[i] for i in eachindex(γ_c) if γ_c[i] == 2.0]

#F2_c = [F_c[i] for i in eachindex(γ_c) if γ_c[i] == X]
#Nu2_c = [Nu_c[i] for i in eachindex(γ_c) if γ_c[i] == X]

# Plot
fig_c = Figure(size=(800, 800))
ax_c = Axis(fig_c[1, 1],
          xlabel="Wind forcing (constant u bc)",
          ylabel="Nu",
          title="Nu vs constant Wind Forcing at Different R/R_c",
          xticks=0:0.001:0.01)

scatter!(ax_c, F1_c, Nu1_c, label="R/R_c = 2")
#scatter!(ax_c, F2_c, Nu2_c, label="R/R_c = 2.5")

axislegend(ax_c)

save("nu_vs_constant_wind_forcing.png", fig_c)

using CairoMakie, CSV, DataFrames

# Load data
df_3 = CSV.read("data/gamma_and_Nu_1.csv", DataFrame)

# Manually compute safe x-limits (for log scale)
xmin = minimum(df_3.gamma)
xmax = maximum(df_3.gamma)

# Create figure and axis with log-log scale
fig_3 = Figure(size=(800, 800))
ax_3 = Axis(fig_3[1, 1];
    xlabel = "γ = R/R_c",
    ylabel = "Nu",
    title = "γ vs Nu (Oceananigans vs Veronis 1966)",
    limits = ((xmin * 0.9, xmax * 1.1), nothing),
    xscale = log10,
    yscale = log10
)

# Plot each source with different style
for src in unique(df_3.source)
    df_src = filter(row -> row.source == src, df_3)
    scatter!(ax_3, df_src.gamma, df_src.Nu; label=src)
end

# Add legend and save
axislegend(ax_3)
save("γ_vs_Nu_comparison.png", fig_3)

#1 - 5-6days low res
#2 - 0-1days high res
#3 - 0-5days high res
#4 - 1-2days very high res -> timed out need to redo at some point. Will take min 24hrs tho

df_4 = CSV.read("data/Changing_Prandtl_number.csv", DataFrame)

# Create figure
fig_4 = Figure(size=(800, 800))
ax_4 = Axis(fig_4[1, 1];
    xlabel = "Prandtl number",
    ylabel = "Nu",
    title = "Nu vs Prandtl number (log-log with best-fit lines)",
    xscale = log10,
    yscale = log10
)

for χ_val_4 in unique(df_4.chi)
    df_subset_4 = filter(row -> row.chi == χ_val_4, df_4)
    scatter!(ax_4, df_subset_4.Prandtl, df_subset_4.Nu; label = "χ = $χ_val_4")

    # Prepare log-transformed DataFrame
    logdf = DataFrame(
        logPr = log10.(df_subset_4.Prandtl),
        logNu = log10.(df_subset_4.Nu)
    )

    # Fit line: logNu ~ logPr
    model = lm(@formula(logNu ~ logPr), logdf)
    intercept = coef(model)[1]
    slope = coef(model)[2]

    # Back-transform to plot in original log-log space
    Pr_range = range(minimum(df_subset_4.Prandtl), stop=maximum(df_subset_4.Prandtl), length=100)
    Nu_fit = 10.0^intercept .* Pr_range .^ slope

    # Plot the fitted line
    lines!(ax_4, Pr_range, Nu_fit; linestyle = :dash, label = "fit χ = $χ_val_4")

    println("χ = $χ_val_4 → Nu ≈ $(round(10^intercept, sigdigits=3)) × Pr^$(round(slope, sigdigits=3))")
end

axislegend(ax_4)

save("Nu_vs_Prandtl_number.png", fig_4)