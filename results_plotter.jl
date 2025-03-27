using CairoMakie

# Data
R = [2, 2, 2, 2, 2, 2, 2, 2, 2,
     2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
F = [1e-10, 0, 1e-8, 1e-7, 5e-7, 2e-7, 1.5e-7, 2.5e-7, 3e-7,
     1e-10, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 5e-7, 2e-7]
Nu = [1.226157172, 1.226157229, 1.225585405, 1.167617251, 1.0, 1.000001444, 1.0862497, 1.0, 1.0,
      1.538779181, 1.538779141, 1.538775225, 1.538383478, 1.497647982, 1.0, 1.0, 1.358498654]

# Separate by R/R_c
F1 = [F[i] for i in eachindex(R) if R[i] == 2.0]
Nu1 = [Nu[i] for i in eachindex(R) if R[i] == 2.0]

F2 = [F[i] for i in eachindex(R) if R[i] == 2.5]
Nu2 = [Nu[i] for i in eachindex(R) if R[i] == 2.5]

# Plot
fig = Figure(resolution=(800, 800))
ax = Axis(fig[1, 1],
          xlabel="Wind forcing (flux)",
          ylabel="Nu",
          title="Nu vs flux Wind Forcing at Different R/R_c")

scatter!(ax, F1, Nu1, label="R/R_c = 2")
scatter!(ax, F2, Nu2, label="R/R_c = 2.5")

axislegend(ax)

# Save figure
save("nu_vs_flux_wind_forcing.png", fig)

#Constant wind forcing
γ_c = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
F_c = [0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 5e-4, 8e-4, 1e-3, 1e-2, 5e-3]
Nu_c = [
    1.7942710918869882,
    1.7942710918869826,
    1.7942710918868796,
    1.7942710918801934,
    1.794271091214404,
    1.7942710246276579,
    1.7942643657080029,
    1.7935959928249512,
    1.775863590460335,
    1.7405832443707834,
    1.7006242785703853,
    1.0,
    1.0
]

F1_c = [F_c[i] for i in eachindex(γ_c) if γ_c[i] == 2.0]
Nu1_c = [Nu_c[i] for i in eachindex(γ_c) if γ_c[i] == 2.0]

#F2_c = [F[i] for i in eachindex(R) if R[i] == 2.5]
#Nu2_c = [Nu[i] for i in eachindex(R) if R[i] == 2.5]

# Plot
fig_c = Figure(resolution=(800, 800))
ax_c = Axis(fig_c[1, 1],
          xlabel="Wind forcing (flux)",
          ylabel="Nu",
          title="Nu vs constant Wind Forcing at Different R/R_c")

scatter!(ax_c, F1_c, Nu1_c, label="R/R_c = 2")
#scatter!(ax_c, F2_c, Nu2_c, label="R/R_c = 2.5")

axislegend(ax_c)

save("nu_vs_constant_wind_forcing.png", fig_c)