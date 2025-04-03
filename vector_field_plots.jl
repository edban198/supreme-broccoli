using CairoMakie

# Parameters
A₀ = 1
α₀ = 1
x₀ = π
y₀ = 0.5

# Vector field functions
u₀(x,y) = A₀ * 2 * (y - y₀) * α₀ * exp(-α₀ * ((x - x₀)^2 + (y - y₀)^2)) -
          A₀ * 2 * (y + y₀) * α₀ * exp(-α₀ * ((x - x₀)^2 + (y + y₀)^2))
v₀(x,y) = - (A₀ * 2 * (x - x₀) * α₀ * exp(-α₀ * ((x - x₀)^2 + (y - y₀)^2)) -
             A₀ * 2 * (x - x₀) * α₀ * exp(-α₀ * ((x - x₀)^2 + (y + y₀)^2)))

# Define a finer grid
x = LinRange(0, 2π, 100)
y = LinRange(-10, 10, 100)
X, Y = [xi for xi in x, _ in y], [yi for _ in x, yi in y]

# Evaluate vector field over the grid
U_vals = u₀.(X, Y)
V_vals = v₀.(X, Y)

# Create figure
fig = Figure(resolution=(1000, 600))
ax = Axis(fig[1, 1]; limits = ((0, 2π), (-10, 10)),
          xlabel = "x", ylabel = "y", title = "Dipole Velocity Field (Streamlines)")

# Plot streamlines
streamplot!(ax, x, y, U_vals, V_vals;
            colormap = :viridis,
            line_width = 2.0,
            arrow_size = 6)

# Save or display
save("OUPUTS/dipole_vec_field.png", fig)