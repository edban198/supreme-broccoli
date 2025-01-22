using Oceananigans
using Oceananigans.Units: minutes, hours, meters, kilometers
using CairoMakie

# Constants
const Nx = 64     # Number of points in horizontal directions
const Nz = 16     # Number of points in the vertical direction
const Lx = 5000     # Domain horizontal extent
const Lz = 1000      # Domain depth
const refinement = 1.2     # Spacing near the surface
const stretching = 12      # Stretching rate at the bottom

# Heaviside function
heaviside(x) = x < 0 ? 0.0 : 1.0

# Convert `Lz` to a Float64
const H = Lz  # Domain depth in meters as Float64

# Sponge region boundaries
sponge_one = -H / 4
sponge_zero = sponge_one + H / 10

# Bottom mask function
function bottom_mask_func(z)
    sponge_one = -H / 4
    sponge_zero = sponge_one + H / 10
    heaviside(-(z - sponge_zero)) * (z - sponge_zero)^2 / (sponge_one - sponge_zero)^2
end

# Generate z-values for plotting
z_values = range(-H, 0, length=500)  # From -H to 0
mask_values = [bottom_mask_func(z) for z in z_values]

# Plot using CairoMakie
fig = Figure(size = (800, 500))
ax = Axis(fig[1, 1],
    title = "Bottom Mask Function",
    xlabel = "z (meters)",
    ylabel = "Mask Value"
)

lines!(ax, z_values, mask_values, label = "Bottom Mask Function")
vlines!(ax, [sponge_one, sponge_zero], color = [:red, :green], linestyle = [:dash, :dash], label = ["sponge_one", "sponge_zero"])
axislegend(ax, position = :rt)
fig