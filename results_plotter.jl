using CairoMakie
using CSV, DataFrames

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