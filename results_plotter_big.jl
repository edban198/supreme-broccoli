#720171 - 256by128withwind
#XXX - 256by128nowind

using CairoMakie
using LaTeXStrings
using CSV
using DataFrames
using GLM
using Statistics

# === Load and prepare data ===
df = CSV.read("data/low_res_big_sim.csv", DataFrame)
df.Prandtl = parse.(Float64, string.(df.Prandtl))
df.chi     = parse.(Float64, string.(df.chi))
df.Nu      = parse.(Float64, string.(df.Nu))
df.taux    = parse.(Float64, string.(df.taux))

R_c = 1100.65
df.R = df.chi .* R_c

# Split datasets
df_nowind = filter(:taux => ==(0.0), df)
df_wind   = filter(:taux => !=(0.0), df)

# === Helper: Plot scatter with line of best fit ===
function plot_scatter_fit(df_subset, xcol, xlabel, groupcol, label_suffix, file_prefix)
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1]; xlabel=xlabel, ylabel="Nu",
              title="Nu vs $xlabel for different $(groupcol) ($label_suffix)",
              xscale=log10, yscale=log10)

    for g in sort(unique(df_subset[!, groupcol]))
        df_group = filter(row -> row[groupcol] == g, df_subset)
        x = df_group[!, xcol]
        y = df_group.Nu

        scatter!(ax, x, y; label="$groupcol = $g")

        model = lm(@formula(log10(y) ~ log10(x)), DataFrame(x=x, y=y))
        a, b = coef(model)
        x_range = range(minimum(x), maximum(x), length=100)
        y_fit = 10^a .* x_range .^ b
        lines!(ax, x_range, y_fit; linestyle=:dash)
    end

    axislegend(ax, position=:rb)
    save("$(file_prefix)_vs_Nu_by_$(xcol).png", fig)
end

# === Plotting ===
plot_scatter_fit(df_nowind, :Prandtl, "Prandtl number", :chi, "No Wind", "no_wind_Pr")
plot_scatter_fit(df_wind, :Prandtl, "Prandtl number", :chi, "With Wind", "with_wind_Pr")
plot_scatter_fit(df_nowind, :R, "Rayleigh number", :Prandtl, "No Wind", "no_wind_R")
plot_scatter_fit(df_wind, :R, "Rayleigh number", :Prandtl, "With Wind", "with_wind_R")

# === Regression summaries ===
function regression_summary(df_subset, label)
    df_log = DataFrame(
        logNu = log.(df_subset.Nu),
        logPr = log.(df_subset.Prandtl),
        logR  = log.(df_subset.R)
    )
    model = lm(@formula(logNu ~ logPr + logR), df_log)
    println("\n📊 Regression Summary for $label")
    display(coeftable(model))
    r2_val = r2(model)
    n = nobs(model)
    k = length(coef(model)) - 1
    adj_r2 = 1 - (1 - r2_val) * (n - 1) / (n - k - 1)
    stderr = sqrt(deviance(model) / dof_residual(model))
    println("R²: ", round(r2_val, digits=4))
    println("Adjusted R²: ", round(adj_r2, digits=4))
    println("Residual Std. Error: ", round(stderr, digits=4))
end

regression_summary(df_nowind, "No Wind")
regression_summary(df_wind, "With Wind")