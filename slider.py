#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
import copy

# Function to update the arrays based on the value of p
def update_arrays(p):
    arrays=copy.deepcopy(currentlist)
    for arr in arrays:
        arr[1] = arr[1] * (1 + 0.5 * arr[1] * np.log(p))
    return arrays

local_range=[0.95,1.05]
with open("point1.0.txt", "r") as file:
    currentlist=eval(file.read())

name1="point"+str(local_range[0])+".txt"
name2="point"+str(local_range[1])+".txt"
with open(name1, "r") as file:
    list1 = eval(file.read())
with open(name2, "r") as file:
    list2 = eval(file.read())

# Function to interpolate between two lists of arrays based on q
def interpolate_lists(list1, list2, q):
    r=(q-local_range[0])/(local_range[1]-local_range[0])
    global currentlist
    currentlist=[ (1 - r) * arr1 + r * arr2 for arr1, arr2 in zip(list1, list2) ]
    return currentlist

# Function to update vib slider, don't affect local slider
def update_vib(event=None):
    p = vib_slider.get()
    updated_arrays = update_arrays(p)
    ax.clear()
    for arr in (updated_arrays):
        ax.plot(arr[0],arr[1],color='r')
    ax.set_title('p = '+str(p))
    ax.set_xlim(0.1,0.5)
    ax.set_ylim(1.0,2.4)
    canvas.draw()

# Function to update local slider, refresh affect vib slider
def update_local(event=None):
    vib_slider.set(1.0)
    q = local_slider.get()
    interpolated_list = interpolate_lists(list1, list2, q)
    ax.clear()
    for arr in interpolated_list:
        ax.plot(arr[0], arr[1],color='r')  # Plot the first row as x (Composition) and the second row as y (T)
    ax.set_title(f'q = {q:.2f}')
    ax.set_xlim(0.1,0.5)
    ax.set_ylim(1.0,2.4)
    canvas.draw()

# Create the main tkinter window
root = tk.Tk()
root.title("Slider Control for p")

# Create a matplotlib figure and axis
fig, ax = plt.subplots()
ax.set_xlabel('Composition')
ax.set_ylabel('T')
ax.set_title('p = 1.00')
ax.set_xlim(0.1,0.5)
ax.set_ylim(1.0,2.4)

# Plot the initial arrays
for arr in (currentlist):
    ax.plot(arr[0],arr[1],color='r')

# Embed the matplotlib figure in the tkinter window
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# Create a slider for p
vib_slider = tk.Scale(root, from_=0.8, to=1.2, resolution=0.01,orient=tk.HORIZONTAL, label="vibration r")
vib_slider.set(1.0)  # Initial value of p
vib_slider.pack()
vib_slider.bind("<Motion>", update_vib)  # Update plot on slider change

local_slider = tk.Scale(root, from_=local_range[0], to=local_range[1], resolution=0.01,orient=tk.HORIZONTAL, label="local parameter p")
local_slider.set(1.0)  # Initial value of p
local_slider.pack()
local_slider.bind("<Motion>", update_local)  # Update plot on slider change

# Run the tkinter main loop
root.mainloop()

