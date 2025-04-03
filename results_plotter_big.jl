using CairoMakie
using LaTeXStrings
using CSV
using DataFrames
using GLM
using Statistics
using Distributions

# === Load and prepare data ===
df = CSV.read("data/high_res_big_sim.csv", DataFrame)
df.Prandtl = parse.(Float64, string.(df.Prandtl))
df.chi     = parse.(Float64, string.(df.chi))
df.Nu      = parse.(Float64, string.(df.Nu))
df.taux    = parse.(Float64, string.(df.taux))

R_c = 1100.65
df.R = df.chi .* R_c

# === Global regression: log(Nu) ~ log(Pr) + log(R) ===
df_log = DataFrame(
    logNu = log.(df.Nu),
    logPr = log.(df.Prandtl),
    logR  = log.(df.R)
)
model = lm(@formula(logNu ~ logPr + logR), df_log)
global_a, global_b, global_c = coef(model)

# === Helper: Plot scatter with global model curves ===
function plot_global_fit(df_subset, xcol, xlabel, groupcol, label_suffix, file_prefix)
    fig = Figure(size=(1000, 800))
    ax = Axis(fig[1, 1];
        xlabel = xlabel,
        ylabel = "Nu",
        title = "Nu vs $xlabel for different $(groupcol) ($label_suffix)",
        xscale = log10, yscale = log10
    )

    for g in sort(unique(df_subset[!, groupcol]))
        df_group = filter(row -> row[groupcol] == g, df_subset)
        x = df_group[!, xcol]
        y = df_group.Nu
        scatter!(ax, x, y; label="$groupcol = $g")

        # Global model line at fixed other param
        x_range = range(minimum(x), maximum(x), length=200)
        if xcol == :Prandtl
            R_val = g * R_c
            y_fit = exp(global_a) .* (x_range .^ global_b) .* (R_val ^ global_c)
        elseif xcol == :R
            Pr_val = g
            y_fit = exp(global_a) .* (x_range .^ global_c) .* (Pr_val ^ global_b)
        end
        lines!(ax, x_range, y_fit; linestyle=:dash)
    end

    axislegend(ax, position=:rb)
    save("OUTPUTS/$(file_prefix)_vs_Nu_by_$(xcol).png", fig)
end

# === Plotting ===
plot_global_fit(df, :Prandtl, "Prandtl number", :chi, "No Wind", "no_wind_Pr")
plot_global_fit(df, :R, "Rayleigh number", :Prandtl, "No Wind", "no_wind_R")

# === Regression summary with Chi-square goodness-of-fit ===
function regression_summary_with_chisq(df_subset, label)
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

    # Chi-square goodness-of-fit (Pearson) for continuous data
    Nu_pred = exp.(predict(model))
    Nu_obs = df_subset.Nu

    # Avoid dividing by zero (filter or clip)
    valid = Nu_pred .> 0
    Nu_obs = Nu_obs[valid]
    Nu_pred = Nu_pred[valid]

    chi2_stat = sum(((Nu_obs .- Nu_pred).^2) ./ Nu_pred)
    dof = length(Nu_obs) - length(coef(model))
    p_val = 1 - cdf(Chisq(dof), chi2_stat)

    println("\n🔍 Pearson Chi-Square Test (continuous approximation):")
    println("Chi² statistic: ", round(chi2_stat, digits=4))
    println("Degrees of freedom: ", dof)
    println("p-value: ", round(p_val, digits=4))
    
end

regression_summary_with_chisq(df, "No Wind")