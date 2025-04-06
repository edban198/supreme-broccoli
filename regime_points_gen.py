import numpy as np
import csv

# ============================================================================
# 1) Define the region in terms of inequalities:
#    Here, x = R (Rayleigh number) and y = Pr (Prandtl number).
#    Adjust these conditions to match the exact boundaries from your table.
# ============================================================================
def in_region(x, y):
    # Example boundaries (please update these as needed):
    # For instance, using:
    #   Condition 1: y < 6.7e-5 * x^(2/3)
    #   Condition 2: y > 4.3e-8 * x^(2/3)
    #   Condition 3: y > 3.5e-8 * x^(1.5)
    #   Condition 4: y < 1.0e22 * x^(-2)
    #   Condition 5: y < 2
    cond1 = y < 6.7e-5 * (x ** (2/3))
    cond2 = y > 4.3e-8 * (x ** (2/3))
    cond3 = y > (3.5 * 10 ** (3/2)) * (x ** (-2))
    cond4 = y < 1.0e22 * (x ** (-2))
    cond5 = y < 2
    return cond1 and cond2 and cond3 and cond4 and cond5

# ============================================================================
# 2) Define the bounding box in log-space and the grid density
#    (These are our initial guesses; the actual region will be determined by in_region.)
# ============================================================================
# For x = R, we want log10(R) from 5 to 10, i.e. R in [1e5, 1e10]
X_MIN = 1207.11172
X_MAX = 1.02908e11
X_MIN_LOG10 = np.log10(X_MIN)  # Convert to log-space
X_MAX_LOG10 = np.log10(X_MAX)  # Convert to log-space

# For y = Pr, we want log10(Pr) from -3 to 1, i.e. Pr in [1e-3, 10^1]
Y_MIN = 0.000030628
Y_MAX = 2
Y_MIN_LOG10 = np.log10(Y_MIN)  # Convert to log-space
Y_MAX_LOG10 = np.log10(Y_MAX)  # Convert to log-space

# Number of steps in each dimension; adjust these to control density:
NUM_X = 12
NUM_Y = 12

# ============================================================================
# 3) Create the grid in log-space and filter the points using in_region
# ============================================================================
x_log_vals = np.linspace(X_MIN_LOG10, X_MAX_LOG10, NUM_X)
y_log_vals = np.linspace(Y_MIN_LOG10, Y_MAX_LOG10, NUM_Y)

points_in_region = []

for lx in x_log_vals:
    x_val = 10 ** lx  # Convert from log-space to actual value
    for ly in y_log_vals:
        y_val = 10 ** ly
        if in_region(x_val, y_val):
            points_in_region.append((x_val, y_val))

print(f"Total grid points: {NUM_X*NUM_Y}")
print(f"Points in region: {len(points_in_region)}")

# ============================================================================
# 4) Compute the minimum and maximum R and Pr among the valid points
# ============================================================================
if points_in_region:
    R_values = [pt[0] for pt in points_in_region]
    Pr_values = [pt[1] for pt in points_in_region]
    R_min, R_max = min(R_values), max(R_values)
    Pr_min, Pr_max = min(Pr_values), max(Pr_values)
    print(f"R range in region: {R_min:.2e} to {R_max:.2e}")
    print(f"Pr range in region: {Pr_min:.2e} to {Pr_max:.2e}")
else:
    print("No points found in the region.")

# ============================================================================
# 5) Write the valid points to a CSV file
# ============================================================================
out_filename = "points_in_region_I.csv"
with open(out_filename, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["R", "Pr"])  # header
    for pt in points_in_region:
        writer.writerow(pt)

print(f"Wrote {len(points_in_region)} points to {out_filename}")