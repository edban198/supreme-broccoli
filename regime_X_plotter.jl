using CSV, DataFrames, GLM, Statistics, CairoMakie, LaTeXStrings, Printf

# --- Function to load and parse data ---
function load_data(filename::String)
    df = CSV.read(filename, DataFrame)
    # Ensure columns are Float64 (assuming CSV columns are named "Pr", "R", "Nu", "τx")
    df.Pr = parse.(Float64, string.(df.Pr))
    df.R  = parse.(Float64, string.(df.R))
    df.Nu = parse.(Float64, string.(df.Nu))
    return df
end

# --- Function to perform linear regression in log–log space ---
function perform_regression(df::DataFrame)
    # Create a new DataFrame with natural logs of Nu, Pr, and R
    df_log = DataFrame(
        logNu = log.(df.Nu),
        logPr = log.(df.Pr),
        logR  = log.(df.R)
    )
    # Fit the linear model: log(Nu) ~ log(Pr) + log(R)
    model = lm(@formula(logNu ~ logPr + logR), df_log)
    return model
end

# --- Function to print the regression summary ---
function print_regression_summary(model)
    ct = coeftable(model)
    r2_val = r2(model)
    n = nobs(model)
    k = length(coef(model)) - 1
    adj_r2 = 1 - (1 - r2_val) * (n - 1) / (n - k - 1)
    stderr = sqrt(deviance(model) / dof_residual(model))
    
    println("\n📊 Regression Summary:")
    display(ct)
    @printf("R²: %.4f\n", r2_val)
    @printf("Adjusted R²: %.4f\n", adj_r2)
    @printf("Residual Std. Error: %.4f\n", stderr)
    
    # Extract coefficients
    global_a = coef(model)[1]
    global_b = coef(model)[2]
    global_c = coef(model)[3]
    pre_factor = exp(global_a)
    println("\nDerived power law:")
    println("Nu ≈ $(round(pre_factor, digits=3)) * Pr^$(round(global_b, digits=3)) * R^$(round(global_c, digits=3))")
end

# --- Function to create a 3D scatter plot ---
function scatter_3d_plot(df::DataFrame, output_filename::String)
    # Compute log values for plotting
    logPr = log.(df.Pr)
    logR  = log.(df.R)
    logNu = log.(df.Nu)
    
    fig = Figure(size = (1000, 800))
    ax = Axis3(fig[1,1],
        xlabel = L"\ln(\mathrm{Pr})",
        ylabel = L"\ln(R)",
        zlabel = L"\ln(\mathrm{Nu})",
        title = "3D Scatter of log(Pr), log(R), log(Nu)"
    )
    scatter!(ax, logPr, logR, logNu; markersize=8, color=:blue)
    save(output_filename, fig)
    println("3D scatter plot saved to ", output_filename)
end

# --- Main script ---
filename = "data/regime_I_data.csv"  # adjust path if needed
df = load_data(filename)

# Perform regression and print summary
model = perform_regression(df)
print_regression_summary(model)

# Create a 3D scatter plot (using log-space values)
scatter_3d_plot(df, "OUTPUTS/regime_I_scatter.png")