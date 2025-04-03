using DataFrames, CairoMakie, CSV, GLM, StatsModels, Printf

df = CSV.read("data/Nu_vs_Pr.csv", DataFrame)

# Separate by χ value
df10 = sort(filter(:chi => ==(10.0), df), :Prandtl)
df15 = sort(filter(:chi => ==(15.0), df), :Prandtl)

# Log–log transforms
df10_log = DataFrame(logPr = log10.(df10.Prandtl), logNu = log10.(df10.Nu))
df15_log = DataFrame(logPr = log10.(df15.Prandtl), logNu = log10.(df15.Nu))

# Linear fits
fit10 = lm(@formula(logNu ~ logPr), df10_log)
fit15 = lm(@formula(logNu ~ logPr), df15_log)

slope10, intercept10 = coef(fit10)[2], coef(fit10)[1]
slope15, intercept15 = coef(fit15)[2], coef(fit15)[1]

A10 = 10^intercept10
A15 = 10^intercept15

@info "χ = 10 → Nu ≈ $(round(A10, sigdigits=3)) * Pr^$(round(slope10, digits=3))"
@info "χ = 15 → Nu ≈ $(round(A15, sigdigits=3)) * Pr^$(round(slope15, digits=3))"

# Prediction lines for plotting
Pr_range = range(minimum(df10.Prandtl), stop=maximum(df10.Prandtl), length=200)
Nu_fit10 = A10 .* Pr_range .^ slope10
Nu_fit15 = A15 .* Pr_range .^ slope15

# Plot
fig = Figure(size=(1000, 600))
ax = Axis(fig[1,1],
    xlabel = "Pr",
    ylabel = "Nu",
    xscale = log10,
    yscale = log10,
    title = "Nu vs Pr (log-log) with power law fits"
)

# Data points
scatter!(ax, df10.Prandtl, df10.Nu, label="χ = 10", color=:blue)
scatter!(ax, df15.Prandtl, df15.Nu, label="χ = 15", color=:red)

# Fit lines
lines!(ax, Pr_range, Nu_fit10, color=:blue, linestyle=:dash, label="Fit χ=10")
lines!(ax, Pr_range, Nu_fit15, color=:red, linestyle=:dash, label="Fit χ=15")

# Annotate slopes
text!(ax, 3, 6, text=@sprintf("slope (χ=10) ≈ %.3f", slope10), color=:blue, align=(:left, :top))
text!(ax, 3, 4.5, text=@sprintf("slope (χ=15) ≈ %.3f", slope15), color=:red, align=(:left, :top))

axislegend(ax, position=:rb)
save("OUTPUTS/Nu_vs_Pr.png", fig)