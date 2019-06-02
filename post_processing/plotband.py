import numpy as np
import matplotlib.pyplot as plt
from extract_data import *

filename1="" #path to the file which contains band structure data
data1=read_file(filename1)

fig=plt.figure()
theplot=fig.add_subplot(111)
theplot.hold(True)
for i in range(1,31):
	theplot.plot(data1[:,0],data1[:,i],'-b')
plt.show()
