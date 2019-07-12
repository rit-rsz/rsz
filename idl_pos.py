import matplotlib.pyplot as plt
import numpy as np

x = []
y = []
f_obj = open('idl_data.csv')
for line in f_obj:
    data = line.split(',')
    x.append(data[0])
    y.append(data[1])

plt.scatter(x, y)
plt.show()
