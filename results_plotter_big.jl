using CSV, DataFrames, GLM, Statistics, Printf, LaTeXStrings, CairoMakie

# ----------------------------------------
# Function to load data
# ----------------------------------------
function load_data(filename::String)
    df = CSV.read(filename, DataFrame)
    # Ensure the columns are parsed as Float64
    df.Prandtl = parse.(Float64, string.(df.Pr))
    df.taux    = parse.(Float64, string.(df.taux))
    df.Nu      = parse.(Float64, string.(df.Nu))
    df.R       = parse.(Float64, string.(df.R))
    # filter out the low-R regime
    df = filter(row -> row.R > 3.5e3, df)
    return df
end

# ----------------------------------------
# Print regression summary for each Pr, including R² and Adjusted R²
# ----------------------------------------
function print_pr_regression_summary(model, pr_val)
    ct    = coeftable(model)
    df_ct = DataFrame(ct)   # first 5 cols: term, estimate, stderr, t, p

    println("\n📊 Regression Summary for Pr ≈ $(pr_val):\n")
    for row in eachrow(df_ct)
        term = row[1]
        est  = row[2]
        se   = row[3]
        tval = row[4]
        pval = row[5]
        @printf("%-10s : %8.4f ± %6.4f, t = %5.2f, p = %6.3g\n",
                term, est, se, tval, pval)
    end

    # Derived power‐law: Nu = A * R^β
    a, b = coef(model)
    A     = exp(a)
    se0   = df_ct[1, 3]                      # std-error of intercept is column 3
    lower = exp(a - 1.96*se0)
    upper = exp(a + 1.96*se0)
    println("\n🧮 Fit: Nu ~ $(round(A,digits=3)) × R^$(round(b,digits=3))")
    println("Prefactor A = $(round(A, sigdigits=4)) [95% CI: $(round(lower,sigdigits=3)) – $(round(upper,sigdigits=3))]")

    # R² and Adjusted R²
    r2_val = r2(model)
    n      = nobs(model)
    k      = length(coef(model)) - 1
    adj_r2 = 1 - (1 - r2_val)*(n - 1)/(n - k - 1)
    @printf("🔎 R²: %.4f, Adjusted R²: %.4f\n", r2_val, adj_r2)
end

# ----------------------------------------
# 3D scatter plot of log(Pr), log(R), log(Nu)
# ----------------------------------------
function plot_3d_scatter(df::DataFrame, output_filename::String; azimuth=0.75π, elevation=0.4)
    logPr = log.(df.Prandtl)
    logR  = log.(df.R)
    logNu = log.(df.Nu)
    fig = Figure(size=(1000,800))
    ax = Axis3(fig[1,1],
        xlabel    = L"\ln(\mathrm{Pr})",
        ylabel    = L"\ln(R)",
        zlabel    = L"\ln(\mathrm{Nu})",
        title     = "3D Scatter: ln(Pr), ln(R), ln(Nu)",
        azimuth   = azimuth,
        elevation = elevation
    )
    sc = scatter!(ax, logPr, logR, logNu;
                  markersize   = 14,
                  color        = logNu,
                  colormap     = :viridis,
                  transparency = true)
    Colorbar(fig[1,2], sc, label="log(Nu)", width=15)

    # add best‐fit lines at Pr=1 and Pr=7
    for (pr_val, colorline) in zip([1.0,7.0], [:red,:blue])
        df_pr = filter(r -> isapprox(r.Prandtl, pr_val; atol=1e-8), df)
        if nrow(df_pr)==0 continue end
        df_pr_log = DataFrame(logR=log.(df_pr.R), logNu=log.(df_pr.Nu))
        model_pr  = lm(@formula(logNu~logR), df_pr_log)
        a, b      = coef(model_pr)
        rmin, rmax = extrema(df_pr.R)
        R_grid    = range(rmin, rmax, length=100)
        logNu_pred = a .+ b .* log.(R_grid)
        lines!(ax,
               fill(log(pr_val), length(R_grid)),
               log.(R_grid),
               logNu_pred;
               color     = colorline,
               linewidth = 3,
               label     = "Pr=$(pr_val): Nu ~ R^$(round(b,digits=3))")
    end
    axislegend(ax)
    save(output_filename, fig)
    println("✨ 3D scatter saved to ", output_filename)
end

# ----------------------------------------
# Animate rotating 3D scatter
# ----------------------------------------
function animate_rotation(df::DataFrame, output_video::String)
    logPr = log.(df.Prandtl)
    logR  = log.(df.R)
    logNu = log.(df.Nu)
    fig = Figure(size=(800,600))
    ax = Axis3(fig[1,1], xlabel="log(Pr)", ylabel="log(R)", zlabel="log(Nu)")
    scatter!(ax, logPr, logR, logNu; color=logNu, colormap=:viridis, markersize=10)
    record(fig, output_video, 1:90) do i
        ax.azimuth[] = i * 2π / 90
    end
    println("🎥 Animation saved to ", output_video)
end

# ----------------------------------------
# Plot Nu vs R for a specific Pr, plus regression summary
# ----------------------------------------
function plot_R_vs_Nu_for_Pr(df::DataFrame, specific_Pr::Float64; tolerance::Float64=0.05)
    df_filtered = filter(r -> abs(r.Prandtl - specific_Pr)/specific_Pr < tolerance, df)
    if nrow(df_filtered)==0
        println("No data for Pr≈", specific_Pr)
        return
    end
    # regression in log–log
    df_log    = DataFrame(logR = log.(df_filtered.R), logNu = log.(df_filtered.Nu))
    model_fit = lm(@formula(logNu ~ logR), df_log)

    # print summary
    print_pr_regression_summary(model_fit, specific_Pr)

    # plot
    fig = Figure(size=(800,600))
    ax  = Axis(fig[1,1],
        xlabel = "Rayleigh Number",
        ylabel = "Nusselt Number",
        xscale = log10,
        yscale = log10,
        xlabelsize = 20,
        ylabelsize = 20,
        xticklabelsize = 20,
        yticklabelsize = 20,
        xticksize = 16,
        yticksize = 16
    )
    scatter!(ax, df_filtered.R, df_filtered.Nu; markersize=16, color=:blue)
    coefs   = coef(model_fit)
    R_range = range(minimum(df_filtered.R), stop=maximum(df_filtered.R), length=100)
    Nu_fit  = exp(coefs[1]) .* R_range .^ coefs[2]
    fit_lbl = "Nu ~ $(round(exp(coefs[1]),digits=3)) × Ra^$(round(coefs[2],digits=3))"
    lines!(ax, R_range, Nu_fit; color=:red, linewidth=2, label=fit_lbl)
    axislegend(ax; labelsize=20)
    out = "OUTPUTS/Nu_vs_R_Pr=$(specific_Pr).png"
    save(out, fig)
    println("📊 Plot saved to ", out)
end

# ----------------------------------------
# Main Script
# ----------------------------------------
filename = "data/report_data_1.csv"
df = load_data(filename)

# 3D scatter & animation
#plot_3d_scatter(df, "OUTPUTS/3d_scatter_with_lines.png"; azimuth=0.75π, elevation=0.4)
#animate_rotation(df, "OUTPUTS/rotating_3d_scatter.mp4")

# Per-Pr regressions & plots
plot_R_vs_Nu_for_Pr(df, 1.0; tolerance=0.05)
plot_R_vs_Nu_for_Pr(df, 7.0; tolerance=0.05)