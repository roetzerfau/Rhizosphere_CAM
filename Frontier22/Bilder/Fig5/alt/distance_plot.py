import matplotlib.pyplot as plt
import csv
import numpy as np

files = []

for p in range(1):
	data =   np.loadtxt(open('porosity_table.csv', 'rb'), delimiter=',')
	#t = np.linspace(0,data.shape[0], data.shape[0])
	dist = [20, 40 ,60 ,80, 100 ,120 , 140]
	labels = ["initial state",  "timestep 50", "root full grown", "shrunken root",  "timestep 666" ,"timestep 999" ]
	t = [0, 50, 99, 196,666,999]
	markers = ["o", "v", "^", "<", ">", "8"]
	colors = plt.cm.Greys(np.linspace(0,1,len(t)))
	#x = list(reader)
	#result = numpy.array(x).astype("float")
	for i in range(len(t)):
		#label = "timestep " + str(t[i] +1)
		#plt.plot(dist, data[t[i]], label = labels[i] )
		plt.plot(dist, data[t[i]], label = labels[i], color = colors[i], marker=markers[i], linestyle='solid' )# color='green',linewidth=2, markersize=12)

	plt.ylim(0,1)
	plt.xlim(0, 150)

	plt.xlabel('Distance from root surface [$\mu$m] ')
	plt.ylabel('Porosity [%]')
	plt.legend()
	plt.savefig('porosity_table_34_normuc_196.png')
	plt.show()
