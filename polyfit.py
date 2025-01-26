#!/usr/bin/python

#multifunctional curve fitting code

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def target_function(x, a, b, c):
    return a + b * np.exp(c * x)

def r_squared(actual, predicted):
    ss_total = np.sum((actual - np.mean(actual)) ** 2)  # Total sum of squares
    ss_residual = np.sum((actual - predicted) ** 2)  # Residual sum of squares
    return 1 - (ss_residual / ss_total)

def fit_local():           #fit local parameter using linear function
    data = np.loadtxt('inv.txt')  # Ensure the file has 3 columns: time, x, y
    p = np.linspace(1.00, 1.09, num=10)  # Third dimension (e.g., time)
    x = data[:, 0]  # x coordinates
    T = data[:, 1]  # y coordinates

    # Step 2: Extract the starting and ending points
    start_point = (p[0], x[0], T[0])
    end_point = (p[-1], x[-1], T[-1])

    # Step 3: Fit a straight line for x and y as functions of t
    # For x: x = m_x * t + c_x
    m_x = (end_point[1] - start_point[1]) / (end_point[0] - start_point[0])  # Slope for x
    c_x = start_point[1] - m_x * start_point[0]  # Intercept for x

    # For y: y = m_y * t + c_y
    m_T = (end_point[2] - start_point[2]) / (end_point[0] - start_point[0])  # Slope for y
    c_T = start_point[2] - m_T * start_point[0]  # Intercept for y

    # Generate predicted x and y values
    x_pred = m_x * p + c_x
    T_pred = m_T * p + c_T

    # Step 4: Compute the residuals (errors) for x and y
    x_residuals = x - x_pred  # Errors in x
    T_residuals = T - T_pred  # Errors in y

    # Compute the mean squared error (MSE) for x and y
    r2_x = r_squared(x, x_pred)  # R^2 for x
    r2_T = r_squared(T, T_pred)  # R^2 for y

    # Print the results
    print(f"Equation for x: x = {m_x:.4f}t + {c_x:.4f}")
    print(f"Equation for y: y = {m_T:.4f}t + {c_T:.4f}")
    print(f"R^2 value for x: {r2_x:.4f}")
    print(f"R^2 value for y: {r2_T:.4f}")

    # Step 5: Visualize the errors
    plt.figure(figsize=(12, 6))

    # Plot for x
    plt.subplot(1, 2, 1)
    plt.scatter(p, x, color='blue', label='Original x')
    plt.plot(p, x_pred, color='red', label='Predicted x')
    plt.xlabel('Local Parameter (p)')
    plt.ylabel('x')
    plt.title(f'Original vs Predicted x (R² = {r2_x:.4f})')
    plt.legend()

    # Plot for T
    plt.subplot(1, 2, 2)
    plt.scatter(p, T, color='green', label='Original T')
    plt.plot(p, T_pred, color='orange', label='Predicted T')
    plt.xlabel('Local Parameter (p)')
    plt.ylabel('Temperature (T)')
    plt.title(f'Original vs Predicted T (R² = {r2_T:.4f})')
    plt.legend()

    plt.tight_layout()
    plt.show()

    # Step 8: Plot the original vs predicted data in x, T space
    plt.figure(figsize=(8, 6))
    plt.scatter(x, T, color='blue', label='Original (x, T)')
    plt.scatter(x_pred, T_pred, color='red', label='Predicted (x, T)')
    plt.plot(x_pred, T_pred, color='green')
    plt.xlabel('composition (x)')
    plt.ylabel('Temperature (T)')
    plt.title('Invariant point prediction in x, T Space')
    plt.legend()
    plt.grid(True)
    plt.show()


# Define the target function: a + b * exp(c * x)

def fit_elastic():
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

if __name__ == '__main__':
    print("this is fitting code")
    fit_elastic()