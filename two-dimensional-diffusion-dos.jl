#install dependencies

#using Pkg
#pkg"add Oceananigans, CairoMakie"

using Oceananigans
using Oceananigans.Units: seconds, hours, day, days

#build grid

grid = RectilinearGrid(GPU(); size = (128,128),
                       x = (-0.5, 0.5), y = (-0.5, 0.5),
                       topology = (Periodic, Periodic, Flat)
)   #assign Flat to x and y for 1d problem

const ν=0.1
const κ=0.1

closure = ScalarDiffusivity(ν=ν, κ=κ)  #use ScalarDiffusivity for either molecular or turbulent diffusion

U(x, y, t) = cos(2π * y)

V(x, y, t) = -cos(2π * x)

model = NonhydrostaticModel(; grid, closure, tracers=:T, background_fields = (u=U, v=V,))    #by default NonhydrostaticModel has no flux b.c on all fields

#set initial condition on temp field
const width = 0.1
initial_temperature(x, y) = exp(-(x^2+y^2) / (2width^2))
set!(model, T=initial_temperature)

#Visualising model data
using CairoMakie
set_theme!(Theme(fontsize = 24, linewidth = 3))

fig = Figure()
axis = (xlabel = "x",
        ylabel = "y")
label = "t=0"
heatmap(model.tracers.T; label, axis)


#Running a simulation - Time-scale for diffusion across a grid cell
min_Δx = minimum_xspacing(model.grid)
diffusion_time_scale = min_Δx^2 / model.closure.κ.T

simulation = Simulation(model,
                        Δt = 30seconds,
                        stop_time = 1day
)

#sim will run for 1000 iterations with a time-step that resolves the time-scale at which our temperature field diffusion_time_scale
#run!(simulation)

#log data to .jld2 file as simulation progresses
filename = "two-dimensional-diffusion"
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

using CairoMakie
set_theme!(Theme(fontsize = 24))

fig = Figure(size = (1000,800))

axis_kwargs = (xlabel = "x", ylabel = "y",
               limits = ((-0.5,0.5), (-0.5,0.5)),
               aspect = AxisAspect(1)
)

using LaTeXStrings

ax_T = Axis(fig[1,1]; title = L"Two-dimensional diffusion in a background field of $\Phi = \frac{1}{2\pi} (\sin(2\pi x / L) + \sin(2\pi y / L))$", axis_kwargs...)

n = Observable(1)

T = @lift T_timeseries[$n]

heatmap!(ax_T, T; colormap = :thermometer) #, colorrange = (0,1)

CairoMakie.activate!(type = "png")
save("initial-temperature.png", fig)

#record movie
frames = 1:length(times)
@info "Making an animation..."
record(fig, filename * ".mp4", frames, framerate=4) do i
    n[] = i
end
