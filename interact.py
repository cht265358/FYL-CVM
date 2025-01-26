import matplotlib.pyplot as plt

# Given matrix
data = [[3, 3], [4, 3], [6, 3], [7.4, 3], [2, 2], [3.4, 2], [6.7, 2], [9.1, 2]]

# Extracting x and y coordinates
x = [point[0] * 0.1 for point in data]
y = [point[1] for point in data]

# Creating the scatter plot
plt.scatter(x, y, color='blue')
plt.ylim(0, 4)
plt.xlim(0, 1)

# Adding labels and title
plt.xlabel('composition')
plt.ylabel('T')
plt.title('example plot for invariant point calculation')

# Display the plot
plt.show()