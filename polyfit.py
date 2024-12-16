#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the target function: a + b * exp(c * x)
def target_function(x, a, b, c):
    return a + b * np.exp(c * x)

print(np.exp(1.5))

with open("elastic_c1.2.txt", "r") as file1:
        # Use eval to parse the file content as a Python list with numpy arrays
            file = eval(file1.read())
            y1y2=file[0]
# Example data (x, y)
x = np.linspace(0,15,16,endpoint=True)
print(x)
y = y1y2[0]

# Fit the curve
initial_guess = [0.39, -0.02, 0.1]  # Initial guesses for a, b, c
params, covariance = curve_fit(target_function, x, y, p0=initial_guess)

# Extract the fitted parameters
a, b, c = params
print(f"Fitted parameters: a = {a:.3f}, b = {b:.3f}, c = {c:.3f}")


# Compute R^2
y_fit = target_function(x, a, b, c)
residuals = y - y_fit
ss_res = np.sum(residuals**2)  # Residual sum of squares
ss_tot = np.sum((y - np.mean(y))**2)  # Total sum of squares
r_squared = 1 - (ss_res / ss_tot)

# Generate the fitted curve
x_fit = np.linspace(min(x), max(x), 500)
y_fit = target_function(x_fit, a, b, c)

# Plot the data and the fitted curve
plt.figure(figsize=(8, 6))
plt.scatter(x, y, color='red', label='Original Data')  # Data points
plt.plot(x_fit, y_fit, color='blue', label=f'Fit: {a:.4f} + {b:.4f} * exp({c:.4f} * x)')  # Fitted curve
plt.xlabel('x')
plt.ylabel('y')
plt.title('Exponential Curve Fit')
s="R^2="+str(r_squared)
plt.text(4,0.35,s)
plt.legend()
plt.grid()
plt.show()