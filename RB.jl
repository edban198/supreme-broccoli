@info"Importing librarys"
#using Pkg
#pkg"add Oceananigans, CairoMakie, LaTeXStrings"

using Printf
using CairoMakie
using LaTeXStrings
using Statistics
using Oceananigans
using Oceananigans.Units: second, seconds, minute, minutes, hour, hours, day, days

const χ = 5
const R = 1707.76 * χ
const Pr = 6.8
const κ = 1e-4
const ν = Pr * κ

filename = "./OUTPUTS/RB_animation"

@info"Setting up model"

const Nx = 1024     # number of points in each of horizontal directions
const Nz = 512          # number of points in the vertical direction

const Lx = 8     # (m) domain horizontal extents
const Lz = 4          # (m) domain depth

grid = RectilinearGrid(CPU(); size = (Nx, Nz),
                       x = (0,Lx),
                       z = (0,Lz),
                       topology = (Periodic, Flat, Bounded)
)

# Buoyancy that depends on temperature:
buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(), constant_salinity=0)

#Set values
const g = buoyancy.gravitational_acceleration
const α = buoyancy.equation_of_state.thermal_expansion
const Δ = ν * κ * R / (g * α * Lz^3)

const t_ff = sqrt(Lz / (g * α * Δ))
const t_ff_days = t_ff / (3600 * 24)
@info "Freefall time in days ~ $t_ff_days"
@info "Freefall time in seconds ~ $t_ff"
const time1 = 12hours  # instead of 24 hours

#Bulk formula
const ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
const u₁₀ = 10    # m s⁻¹, average wind velocity 10 meters above the ocean
const cᴰ = 2.5e-3 # dimensionless drag coefficient
const ρₐ = 1.225  # kg m⁻³, average density of air at sea-level
const τx = 0#(κ/Lz)^2 * ρₒ # m² s⁻²

T_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(0), bottom = ValueBoundaryCondition(Δ))

u_bcs = FieldBoundaryConditions(top = ValueBoundaryCondition(0), bottom = ValueBoundaryCondition(0))

closure = ScalarDiffusivity(ν=ν,κ=κ)

const f = 10 * κ / Lz^2

model = NonhydrostaticModel(; grid, buoyancy,
                            tracers = (:T),
                            closure = closure,
                            boundary_conditions = (T=T_bcs,u=u_bcs,)
)

# Initial conditions

# Random noise
Ξ(x,z) = randn()

# Temperature initial condition: a stable density gradient with random noise superposed.
noise_amplitude = 1e-3
Tᵢ(x, z) = Δ * (1 - z/Lz) + noise_amplitude * Ξ(x, z)
uᵢ(x, z) = noise_amplitude * Ξ(x, z)

# set the model fields using functions or constants:
set!(model, u=uᵢ, w=uᵢ, T=Tᵢ)

# Setting up sim

simulation = Simulation(model, Δt=0.005, stop_time=time1)
wizard = TimeStepWizard(cfl=0.5, max_Δt=0.05)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(50))

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|u|) = %.1e ms⁻¹, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.u), maximum(abs, sim.model.velocities.w),
                                prettytime(sim.run_wall_time)
)

add_callback!(simulation, progress_message, IterationInterval(100))

# OUTPUTS
outputs = (
    w = model.velocities.w,
    u = model.velocities.u,
    T = model.tracers.T,
    #avg_T = mean(model.tracers.T, dims=(1,2)),
    #s = sqrt(model.velocities.u^2 + model.velocities.w^2)
)

const data_interval = 5minutes

simulation.output_writers[:full_outputs] = JLD2OutputWriter(
    model, outputs,
    schedule = TimeInterval(data_interval),
    filename = filename * ".jld2",
    overwrite_existing = true
)

run!(simulation)
@info"Plotting animation"

T_timeseries = FieldTimeSeries(filename * ".jld2", "T")
#s_timeseries = FieldTimeSeries(filename * ".jld2", "s")
w_timeseries = FieldTimeSeries(filename * ".jld2", "w")
u_timeseries = FieldTimeSeries(filename * ".jld2", "u")
#avg_T_timeseries = FieldTimeSeries(filename * ".jld2", "avg_T")
times = T_timeseries.times

set_theme!(Theme(fontsize = 24))

fig = Figure(size = (800,1200))

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)",
               aspect = DataAspect()
)

ax_T = Axis(fig[2,1]; title = L"Temperature, $T$", axis_kwargs...)
#ax_s = Axis(fig[3,1]; title = L"Speed, $s = \sqrt{u^2+v^2}$", axis_kwargs...)
ax_w = Axis(fig[3,1]; title = L"Vertical Velocity, $w$", axis_kwargs...)
ax_u = Axis(fig[4,1]; title = L"Horizontal Velocity, $u$", axis_kwargs...)
#ax_avg_T = Axis(fig[2,3]; title = L"Average Temperature over $x$", xlabel = "T", ylabel = "z(m)")

n = Observable(1)

T = @lift T_timeseries[$n]
#s = @lift s_timeseries[$n]
w = @lift w_timeseries[$n]
u = @lift u_timeseries[$n]
#avg_T = @lift vec(dropdims(avg_T_timeseries[:, :, :, $n], dims=(1,2)))

Tlims = (minimum(abs, interior(T_timeseries)), maximum(abs, interior(T_timeseries)))
#slims = (minimum(abs, interior(s_timeseries)), maximum(abs, interior(s_timeseries)))
wlims = (minimum(interior(w_timeseries)), maximum(abs, interior(w_timeseries)))
ulims = (minimum(interior(u_timeseries)), maximum(abs, interior(u_timeseries)))
#=
xlims!(ax_avg_T, Tlims)
ylims!(ax_avg_T, 0, Lz)
z_vec = LinRange(0, Lz, Nz)
lines!(ax_avg_Tavg_T, color=:red)
=#
hm_T = heatmap!(ax_T, T; colormap = :thermometer, colorrange = Tlims)
Colorbar(fig[2,2], hm_T)
#hm_s = heatmap!(ax_s, s; colormap = :speed, colorrange = slims)
#Colorbar(fig[3,2], hm_s)
hm_w = heatmap!(ax_w, w; colormap = :speed, colorrange = wlims)
Colorbar(fig[3,2], hm_w)
hm_u = heatmap!(ax_u, u; colormap = :speed, colorrange = ulims)
Colorbar(fig[4,2], hm_u)
#=
using Interpolations

zrange = axes(w_timeseries, 3)  # the actual valid indices in the 3rd dimension
w_center_timeseries = 0.5 .* (
    w_timeseries[:, :, zrange[1:end-1], :] .+ w_timeseries[:, :, zrange[2:end], :]
)

wT_timeseries = w_center_timeseries .* T_timeseries

wT_avg_timeseries = mean(wT_timeseries, dims=(1))

wT_avg = @lift vec(dropdims(wT_avg_timeseries[:, :, :, $n], dims=(1,2)))

ax_wT = Axis(fig[3,3];
    title = L"Vertical Heat Flux, $wT$ (averaged over $x$)",
    xlabel = "wT", ylabel = "z (m)"
)

wTlims = (minimum(wT_avg_timeseries), maximum(wT_avg_timeseries))
xlims!(ax_wT, wTlims)
ylims!(ax_wT, 0, Lz)

lines!(ax_wT, wT_avg, color=:red)
=#
@info "calculating Nusselt number"
#Nusselt num
w_timeseries = FieldTimeSeries(filename * ".jld2", "w")
zrange = axes(w_timeseries, 3)  # the actual valid indices in the 3rd dimension
w_center_timeseries = 0.5 .* (
    w_timeseries[:, :, zrange[1:end-1], :] .+ w_timeseries[:, :, zrange[2:end], :]
)

wT_timeseries = w_center_timeseries .* T_timeseries

avg_wT = mean(wT_timeseries)

Nu = 1 + (Lz / (κ * Δ)) * avg_wT

@info "τx = $τx"
@info "Pr = $Pr"
@info "Nu = $Nu"
@info "R = $R"
@info "R/R_c = $χ"
@info "data for csv: $Pr,$R,$Nu,$τx"

title = @lift "t = " * prettytime(times[$n])
Label(fig[1, :], title, fontsize = 24, tellwidth=true)

#record movie
frames = 1:length(times)
@info "Making an animation..."
record(fig, filename * ".mp4", frames, framerate=16) do i
    n[] = i
end

function save_snapshot_at_time(desired_time, output_filename::String="snapshot.png")
    # find the index of the time closest to desired_time
    _, idx = findmin(abs.(times .- desired_time))
    # update the Observable so the figure redraws that frame
    n[] = idx
    # save the current fig to PNG
    save(output_filename, fig)
    @info "Saved snapshot at time=$desired_time (index=$idx) to $output_filename"
end

# Example usage: save the frame closest to t = 6 hours
save_snapshot_at_time(time1, "OUTPUTS/RB_snapshot.png")

if isfile(filename * ".jld2")
    rm(filename * ".jld2"; force=true)
    @info "Deleted file: " * (filename * ".jld2")
end