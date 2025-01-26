#!/usr/bin/python

import numpy as np
import os
import freeenergy

# Example usage
prob = np.random.rand(2, 2, 2, 2)  # Replace with your actual probabilities
energy = freeenergy.compute_binary_cluster_energy(prob, lattice="BCC")
print("Computed Energy:", energy)