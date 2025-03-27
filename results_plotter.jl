using Plots

# Data
R = [2, 2, 2, 2, 2, 2, 2, 2, 2,
     2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
F = [1e-10, 0, 1e-8, 1e-7, 5e-7, 2e-7, 1.5e-7, 2.5e-7, 3e-7,
     1e-10, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 5e-7, 2e-7]
Nu = [1.226157172, 1.226157229, 1.225585405, 1.167617251, 1.0, 1.000001444, 1.0862497, 1.0, 1.0,
      1.538779181, 1.538779141, 1.538775225, 1.538383478, 1.497647982, 1.0, 1.0, 1.358498654]

# Separate groups
F1 = F[[i for i in eachindex(R) if R[i] == 2.0]]
Nu1 = Nu[[i for i in eachindex(R) if R[i] == 2.0]]

F2 = F[[i for i in eachindex(R) if R[i] == 2.5]]
Nu2 = Nu[[i for i in eachindex(R) if R[i] == 2.5]]

# Plot
scatter(F1, Nu1, label="R/R_c = 2", xscale=:log10, xlabel="Wind forcing (flux)", ylabel="Nu", legend=:topright)
scatter!(F2, Nu2, label="R/R_c = 2.5")
