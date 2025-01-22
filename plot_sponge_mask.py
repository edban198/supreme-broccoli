import numpy as np
import matplotlib.pyplot as plt

# Constants
H = 1000  # Domain depth in meters

# Heaviside function
def heaviside(x):
    return np.where(x < 0, 0.0, 1.0)

# Sponge region boundaries
sponge_one = -H / 2
sponge_zero = sponge_one + H / 10

# Bottom mask function
def bottom_mask_func(z):
    return heaviside(-(z - sponge_zero)) * 200 * ((z - sponge_zero) ** 2) / ((sponge_one - sponge_zero) ** 2)

# Generate z-values for plotting
z_values = np.linspace(-H, 0, 500)  # From -H to 0 with 500 points
mask_values = bottom_mask_func(z_values)

# Plot the bottom mask function with flipped axes
plt.figure(figsize=(10, 6))
plt.plot(mask_values, z_values, label="Bottom Mask Function")  # Flip axes
plt.axhline(sponge_one, color='red', linestyle='--', label="sponge_one")
plt.axhline(sponge_zero, color='green', linestyle='--', label="sponge_zero")
plt.title("Bottom Mask Function")
plt.xlabel("Mask Value")
plt.ylabel("z (meters)")
plt.legend()
plt.grid()

# Save the plot as an image
plt.savefig("OUTPUTS/bottom_mask_function_flipped.png", dpi=300, bbox_inches='tight')
plt.show()

print("Plot saved as 'OUTPUTS/bottom_mask_function_flipped.png'")