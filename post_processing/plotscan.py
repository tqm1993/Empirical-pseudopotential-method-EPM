from numpy import *
from extract_data import *
import matplotlib.pyplot as plt


filename="" #path to the file which contains energy scan data
data=read_file(filename)


col_tick=[6,8]

col_ls=['.r','.r','.b','.b','-r','-r','-b','-b']
fig=plt.figure()
theplot=fig.add_subplot(111)
theplot.hold(True)

for i in col_tick:
	theplot.plot(data[:,0],data[:,i],col_ls[i-1])
	
ax=plt.gca()

plt.show()
