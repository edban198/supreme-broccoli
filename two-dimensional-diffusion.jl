#install dependencies

#using Pkg
#pkg"add Oceananigans, CairoMakie"

using CairoMakie
using LaTeXStrings
using Oceananigans
using Oceananigans.Units: meter, meters, day, days, hours, hour, minutes, minute, seconds
using Printf

#build grid

filename = "OUTPUTS/2d_diffusion"

const Lx = 1
const Ly = 1
const Nx = 512
const Ny = 512

grid = RectilinearGrid(GPU(); size=(Nx, Ny),
                       x=(-Lx/2, Lx/2), y=(-Ly/2, Ly/2),
                       topology=(Periodic, Periodic, Flat)
)
   #assign Flat to x and y for 1d problem

const ν=0
const κ=1e-6

closure = ScalarDiffusivity(ν=ν, κ=κ)  #use ScalarDiffusivity for either molecular or turbulent diffusion

const A = 1

model = NonhydrostaticModel(; grid,
                            closure,
                            advection = UpwindBiasedFifthOrder(),
                            tracers=:T
)    #by default NonhydrostaticModel has no flux b.c on all fields

#set initial condition on temp field
const width = 0.1
const A₀ = 1
initial_temperature(x, y) = A₀ * exp(-(x^2+y^2) / (2width^2))
# Use correct grid locations for velocities
U(x, y) = A * cos(2π * y)  # u-velocity at x-face centers
V(x, y) = -A * cos(2π * x) # v-velocity at y-face centers

# Initialize velocities with correct field types
uᵢ = XFaceField(grid)
vᵢ = YFaceField(grid)

set!(uᵢ, U)
set!(vᵢ, V)

set!(model, u=uᵢ, v=vᵢ, T=initial_temperature)

#Visualising model data
using CairoMakie
set_theme!(Theme(fontsize = 24, linewidth = 3))

simulation = Simulation(model,
                        Δt = 0.005seconds,
                        stop_time = 8seconds
)

wizard = TimeStepWizard(cfl=0.4, max_change=1.1, max_Δt=0.001seconds)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|u|) = %.1e ms⁻¹, max(|v|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.u), maximum(abs, sim.model.velocities.v),
                                prettytime(sim.run_wall_time)
)
add_callback!(simulation, progress_message, IterationInterval(20))

#log data to .jld2 file as simulation progresses
simulation.output_writers[:temperature] = JLD2OutputWriter(model, model.tracers,
                                                           filename = filename * ".jld2",
                                                           schedule = TimeInterval(0.05seconds),
                                                           overwrite_existing = true
)

run!(simulation)

#animate results

T_timeseries = FieldTimeSeries(filename * ".jld2", "T")
times = T_timeseries.times

# Convert everything to regular arrays upfront
x = Array(xnodes(grid, Center()))
y = Array(ynodes(grid, Center()))
T_data = Array(interior(T_timeseries))
Tlims = extrema(T_data)

set_theme!(Theme(fontsize = 24, linewidth = 3))

fig = Figure(size = (1000, 800), fontsize=24)

ax = Axis(fig[1, 1],
    xlabel = L"x \; (m)",
    ylabel = L"y \; (m)",
    aspect = DataAspect(),
    limits = ((-Lx/2, Lx/2), (-Ly/2, Ly/2)),
    xgridvisible = false,
    ygridvisible = false,
    xlabelsize = 24,
    ylabelsize = 24,
    xticklabelsize = 24,
    yticklabelsize = 24,
    xticksize = 16,
    yticksize = 16,
    xticks = -0.5:0.5:0.5,
    yticks = -0.5:0.5:0.5
)

n = Observable(1)
T = @lift T_data[:, :, 1, $n]

Tlims = (minimum(interior(T_timeseries)), maximum(interior(T_timeseries)))

xT, yT = nodes(T_timeseries)

hm = heatmap!(ax, x, y, T;
    colormap=:thermometer,
    colorrange=Tlims,
    interpolate=true  # For smoother visualization
)

tickvals  = range(Tlims[1], Tlims[2], length = 6)
ticklabels = string.(round.(tickvals, digits = 2))

Colorbar(fig[1, 2], hm,
    label = L"T",
    width = 12,
    labelsize = 24,
    ticks = (tickvals, ticklabels)
)

title = @lift "t = " * string(round(times[$n], digits=2))

Label(fig[0, :], title, fontsize=24, tellwidth=false)

frames = 1:length(times)

@info "Making a neat animation of vorticity and speed..."

record(fig, filename * ".mp4", frames, framerate=8) do i
    n[] = i
end