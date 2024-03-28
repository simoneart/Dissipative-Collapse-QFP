import numpy as np
import matplotlib.pyplot as plt

# Python code to read data saved by Fortran program

# Open file for reading
with open('solt.txt', 'r') as f:
    lines = f.readlines()

# Initialize arrays
array1 = []
array2 = []

# Parse lines to extract data
for line in lines:
    values = line.split()
    array1.append(float(values[0]))
    array2.append(float(values[1]))

# Now array1 and array2 contain your data

plt.figure(1)
plt.grid()
plt.plot(array1, array2)
plt.ylabel('f')
plt.xlabel('x')
plt.title("Snapshot of the solution")
