import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
import time

def update_plot(val):
    energy = float(scale.get())
    elastic = float(scale2.get())
    vib = float(scale3.get())
    y = ysave+(vib-1)*ysave**2
    #y = energy * ysave
    scatter.set_offsets(np.c_[x, y])  # Update x and y data
    ax.relim()
    ax.autoscale_view()
    canvas.draw()

# Data
x = np.array([])
y = np.array([])
with open("10.txt", "r") as file:
    # Use eval to parse the file content as a Python list with numpy arrays
    xt = eval(file.read())
    for i in range(len(xt)):
        x=np.append(x,xt[i][0])
        x=np.append(x,xt[i][1])
        y=np.append(y,xt[i][2])
        y=np.append(y,xt[i][2])
    ysave=np.copy(y)

# Tkinter setup
root = tk.Tk()
root.title("Interactive Plot")

# Matplotlib figure
fig, ax = plt.subplots()

scatter = ax.scatter(x, y)
ax.set_title("Interactive Plot")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_ylim(0.0,3.0)

canvas = FigureCanvasTkAgg(fig, master=root)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack()

# Slider
scale = tk.Scale(root, from_=0.5, to=2.0, resolution=0.1, orient=tk.HORIZONTAL, label="energy")
scale.pack()
scale.set(1)
scale.bind("<Motion>", update_plot)  # Update plot on slider change

scale2 = tk.Scale(root, from_=-3.0, to=3.0, resolution=0.1, orient=tk.HORIZONTAL, label="elastic")
scale2.pack()
scale2.set(0)
scale2.bind("<Motion>", update_plot)  # Update plot on slider change

scale3 = tk.Scale(root, from_=0.5, to=1.5, resolution=0.1, orient=tk.HORIZONTAL, label="vibration")
scale3.pack()
scale3.set(1)
scale3.bind("<Motion>", update_plot)  # Update plot on slider change

root.mainloop()