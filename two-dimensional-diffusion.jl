#install dependencies

#using Pkg
#pkg"add Oceananigans, CairoMakie"

using CairoMakie
using LaTeXStrings
using Oceananigans
using Oceananigans.Units: meter, meters, day, days, hours, hour, minutes, minute, seconds
using Printf

#build grid

filename = "2d_diffusion"

const Lx = 10
const Ly = 10
const Nx = 64
const Ny = 64

grid = RectilinearGrid(CPU(); size = (Nx,Ny),
                       x = (-Lx/2, Lx/2), y = (-Ly/2, Ly/2),
                       topology = (Periodic, Periodic, Flat)
)   #assign Flat to x and y for 1d problem

const ν=0.1
const κ=0.1

closure = ScalarDiffusivity(ν=ν, κ=κ)  #use ScalarDiffusivity for either molecular or turbulent diffusion

const A = 2π 

U(x, y, t) = A * cos(2π * y)

V(x, y, t) = -A * cos(2π * x)

model = NonhydrostaticModel(; grid,
                            closure,
                            tracers=:T,
                            background_fields = (u=U, v=V,)
)    #by default NonhydrostaticModel has no flux b.c on all fields

#set initial condition on temp field
const width = 0.5
const A₀ = 5
initial_temperature(x, y) = A₀ * exp(-(x^2+y^2) / (2width^2))
set!(model, T=initial_temperature)

#Visualising model data
using CairoMakie
set_theme!(Theme(fontsize = 24, linewidth = 3))

simulation = Simulation(model,
                        Δt = 30seconds,
                        stop_time = 1day
)

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|u|) = %.1e ms⁻¹, max(|v|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.u), maximum(abs, sim.model.velocities.v),
                                prettytime(sim.run_wall_time)
)

add_callback!(simulation, progress_message, IterationInterval(100))

#log data to .jld2 file as simulation progresses
simulation.output_writers[:temperature] = JLD2OutputWriter(model, model.tracers,
                                                           filename = filename * ".jld2",
                                                           schedule = IterationInterval(10minutes),
                                                           overwrite_existing = true
)

run!(simulation)

#animate results

T_timeseries = FieldTimeSeries(filename * ".jld2", "T")
times = T_timeseries.times

xT, yT = nodes(T_timeseries)

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (1000,800))

axis_kwargs = (xlabel = "x", ylabel = "y",
               limits = ((-Lx/2, Lx/2), (-Ly/2, Ly/2)),
               aspect = AxisAspect(1)
)


ax_T = Axis(fig[1,1]; title = L"Two-dimensional diffusion in a background field of $\Phi = \frac{1}{2\pi} (\sin(2\pi x / L) + \sin(2\pi y / L))$", axis_kwargs...)

n = Observable(1)

T = @lift T_timeseries[$n]

Tlims = (minimum(abs, interior(T_timeseries)), maximum(abs, interior(T_timeseries)))

hm_T = heatmap!(ax_T, T; colormap = :thermometer, colorrange = Tlims)
Colorbar(fig[1,2], hm_T)

title = @lift "t = " * prettytime(times[$n])
Label(fig[1, 1:2], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
record(fig, filename * ".mp4", frames, framerate=4) do i
    n[] = i
end