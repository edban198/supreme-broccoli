using DataFrames, CairoMakie, CSV, GLM, StatsModels, Printf, LaTeXStrings, Colors

# ---------------- Helper Function for Log–Scale Limits ----------------
function log_margin_limits(x::AbstractVector{<:Real}; margin_factor=1.2)
    # Filter out nonpositive values (required for log scale)
    x_pos = filter(v -> v > 0, x)
    if isempty(x_pos)
        return (1e-3, 1.0)
    end
    local_min = minimum(x_pos)
    local_max = maximum(x_pos)
    return (local_min / margin_factor, local_max * margin_factor)
end

# ---------------- Data Loading and Cleaning ----------------
function load_and_clean_data(file_path::String)
    df = CSV.read(file_path, DataFrame)
    df = dropmissing(df)
    # Keep only rows where Pr > 0 (required for log-scale)
    return filter(row -> row.Pr > 0, df)
end

# ---------------- Return Sorted Unique Rayleigh Values ----------------
function get_unique_rayleigh(df::DataFrame)
    return sort(unique(df.R))
end

# ---------------- Regression and Plotting Function ----------------
"""
    fit_and_plot_by_R!(ax, df, R_values; palette, filter_nu=true)

For each Rayleigh number in `R_values`:
 - Plots all data points from `df`
 - Computes a best-fit line in log–log space using
   either all points (if `filter_nu==false`) or only points where Nu > 1 (if true)

Returns a DataFrame with fitted parameters (R, A, and β).
"""
function fit_and_plot_by_R!(ax, df::DataFrame, R_values; palette, filter_nu::Bool=true)
    results = DataFrame(R=Float64[], A=Float64[], β=Float64[])
    for (i, R) in enumerate(R_values)
        # Filter data for this Rayleigh number and sort by Pr.
        df_R = filter(row -> row.R == R, df)
        df_R = sort(df_R, :Pr)
        
        # Plot all data points (always plot them, regardless of Nu value)
        scatter!(ax, df_R.Pr, df_R.Nu,
                 label = "R = $(round(R, digits=0))",
                 color = palette[i])
        
        # For regression, choose filtering based on the flag:
        if filter_nu
            df_R_fit = filter(row -> row.Nu > 1, df_R)
        else
            df_R_fit = df_R
        end
        
        # Skip regression if fewer than 2 data points
        if nrow(df_R_fit) < 2
            @warn "Skipping regression for R = $R (fewer than 2 points after filtering)."
            continue
        end
        
        # Fit a model in log–log space: log(Nu) = log(A) + β·log(Pr)
        logPr = log10.(df_R_fit.Pr)
        logNu = log10.(df_R_fit.Nu)
        df_R_log = DataFrame(logPr=logPr, logNu=logNu)
        model = lm(@formula(logNu ~ logPr), df_R_log)
        intercept, slope = coef(model)[1], coef(model)[2]
        A = 10.0^intercept
        
        # Create a LaTeX annotation string.
        eqn_str = @sprintf("R = %.0f: Nu ≈ %.3f · Pr^{%.3f}", R, A, slope)
        eqn = LaTeXString(eqn_str)
        
        # Generate the regression line over the range of Pr in the fit data.
        Pr_range = range(minimum(df_R_fit.Pr), stop=maximum(df_R_fit.Pr), length=200)
        Nu_fit = A .* (Pr_range .^ slope)
        lines!(ax, Pr_range, Nu_fit, color=palette[i], linestyle=:dash)
        
        # Annotate the axis with the fitted equation.
        y_pos = 10^(log10(maximum(df_R.Nu)) - 0.2 * i)
        text!(ax, minimum(df_R.Pr)*1.05, y_pos,
              text = eqn, color=palette[i], align=(:left, :center))
        
        @info "R = $R → Nu ≈ $(round(A, sigdigits=3)) · Pr^$(round(slope, digits=3))"
        push!(results, (R=R, A=A, β=slope))
    end
    return results
end

# ---------------- Single-Condition Plot ----------------
function plot_nu_vs_pr_condition(df::DataFrame, condition_title::String; filter_nu::Bool=true)
    # Compute axis limits with a margin.
    (x_min, x_max) = log_margin_limits(df.Pr; margin_factor=2.0)
    (y_min, y_max) = log_margin_limits(df.Nu; margin_factor=2.0)
    
    # Create figure and axis.
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1,1],
              xlabel = "Pr",
              ylabel = "Nu",
              title = "Nu vs Pr " * condition_title,
              xscale = log10, yscale = log10)
    ax.limits[] = ((x_min, x_max), (y_min, y_max))
    
    # Prepare a color palette based on the number of unique R values.
    R_values = get_unique_rayleigh(df)
    palette = [
        :dodgerblue, :crimson, :forestgreen, :darkorange, :mediumvioletred,
        :teal, :indigo, :gold, :darkviolet, :mediumseagreen, :orangered, :slateblue,
        :darkslategray, :saddlebrown, :darkkhaki, :lightcoral, :lightsteelblue,
        :lightpink, :lightgreen, :lightyellow, :lightgray, :lightblue, :lightcyan,
        :lightgoldenrod
    ][1:length(R_values)]
    
    # Call the regression routine with the chosen filtering flag.
    results = fit_and_plot_by_R!(ax, df, R_values; palette=palette, filter_nu=filter_nu)
    axislegend(ax, position = :rb)
    return fig, results
end

# ---------------- Side-by-Side Comparison Plot ----------------
function plot_comparison(data_no::DataFrame, data_with::DataFrame, output_path::String)
    fig = Figure(size = (1600, 600))
    
    # Compute shared y-limits from the union of both datasets.
    allNu = vcat(data_no.Nu, data_with.Nu)
    (y_min_all, y_max_all) = log_margin_limits(allNu; margin_factor=2.0)
    
    # "No Wind" subplot – use all data for regression (filter_nu = false)
    let
        (x_min_no, x_max_no) = log_margin_limits(data_no.Pr; margin_factor=2.0)
        ax_no = Axis(fig[1, 1],
                     xlabel = "Pr", ylabel = "Nu",
                     title = "Nu vs Pr (No Wind)",
                     xscale = log10, yscale = log10)
        ax_no.limits[] = ((x_min_no, x_max_no), (y_min_all, y_max_all))
        Rvals_no = get_unique_rayleigh(data_no)
        palette_no = [:dodgerblue, :crimson, :forestgreen, :darkorange]  # Adjust palette as needed
        # For no wind, use all data points for the best-fit line.
        fit_and_plot_by_R!(ax_no, data_no, Rvals_no; palette=palette_no, filter_nu=false)
        axislegend(ax_no, position = :rb)
    end
    
    # "Wind Stress" subplot – keep the Nu > 1 filter.
    let
        (x_min_w, x_max_w) = log_margin_limits(data_with.Pr; margin_factor=2.0)
        ax_with = Axis(fig[1, 2],
                       xlabel = "Pr", ylabel = "Nu",
                       title = "Nu vs Pr (Wind Stress)",
                       xscale = log10, yscale = log10)
        ax_with.limits[] = ((x_min_w, x_max_w), (y_min_all, y_max_all))
        Rvals_w = get_unique_rayleigh(data_with)
        palette_w = [:dodgerblue, :crimson, :forestgreen, :darkorange]
        # For wind, use the filter (default: filter_nu = true).
        fit_and_plot_by_R!(ax_with, data_with, Rvals_w; palette=palette_w, filter_nu=true)
        axislegend(ax_with, position = :rb)
    end
    
    save(joinpath(output_path, "Nu_vs_Pr_Comparison.png"), fig)
end

# ---------------- Top-Level Code ----------------
data_file = "data/report_data_f2.csv"  # CSV file with columns: Pr, R, Nu, taux, chi
output_path = "OUTPUTS"
mkpath(output_path)

data = load_and_clean_data(data_file)
data_no_wind = filter(row -> row.taux == 0.0, data)
data_with_wind = filter(row -> row.taux == 0.00001603125, data)

# Produce the side-by-side plot with shared y–limits.
plot_comparison(data_no_wind, data_with_wind, output_path)