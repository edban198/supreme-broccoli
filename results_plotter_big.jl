using CairoMakie
using LaTeXStrings
using CSV
using DataFrames
using GLM
using Statistics

# Load the data
df = CSV.read("data/low_res_big_sim.csv", DataFrame)

# Convert columns if needed
df.Prandtl = parse.(Float64, string.(df.Prandtl))
df.chi = parse.(Float64, string.(df.chi))
df.Nu = parse.(Float64, string.(df.Nu))

# Create figure
fig = Figure(size = (1000, 800))
ax = Axis(fig[1, 1];
    xlabel = L"\chi = R / R_c",
    ylabel = "Nu",
    title = "Nu vs χ for different Prandtl numbers",
    xscale = log10,
    yscale = log10,
)

# Loop over unique Pr values and plot χ vs Nu and best-fit line
for Pr in sort(unique(df.Prandtl))
    df_subset = filter(row -> row.Prandtl == Pr, df)
    sort!(df_subset, :chi)

    # Scatter plot
    scatter!(ax, df_subset.chi, df_subset.Nu; label = "Pr = $Pr")

    # Prepare log-transformed data for regression
    log_chi = log10.(df_subset.chi)
    log_Nu = log10.(df_subset.Nu)
    reg_df = DataFrame(logχ = log_chi, logNu = log_Nu)

    # Fit model: logNu ~ logχ
    model = lm(@formula(logNu ~ logχ), reg_df)
    a = coef(model)[1]     # intercept (log10(a))
    b = coef(model)[2]     # slope

    # Generate fitted line
    χ_range = range(minimum(df_subset.chi), maximum(df_subset.chi), length=100)
    Nu_fit = 10^a .* χ_range .^ b

    lines!(ax, χ_range, Nu_fit; linestyle = :dash)
    
    @info "Pr = $Pr → Nu ≈ $(round(10^a, sigdigits=3)) × χ^$(round(b, sigdigits=3))"
end

axislegend(ax, position = :rb)
save("χ_vs_Nu_by_Pr.png", fig)

# BIG model regression
df.R = 1100.65 .* df.chi
df.logR = log10.(df.R)
df.logPr = log10.(df.Prandtl)
df.logNu = log10.(df.Nu)

model = lm(@formula(logNu ~ logR + logPr), df)

println("regression summary:")
display(coeftable(model))