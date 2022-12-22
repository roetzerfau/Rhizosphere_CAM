import matplotlib.pyplot as plt
import csv
import numpy as np

files = ['agg_19','agg_19_mucilage','agg_34','agg_34_mucilage' ]
names = ['a) 18 % clay', 'b) 18 % clay + mucilage', 'c) 33 % clay', 'd) 33 % clay + mucilage']
dist = [20, 40 ,60 ,80, 100 ,120 , 140]
labels = ["initial state",  "75 days", "100 days", "150 days",  "500 days" ,"1000 days" ]
t = [0, 75, 100,150, 500, 999]
fig, axs = plt.subplots(2,2, figsize=(20, 10))
fig.subplots_adjust(hspace = .25, wspace=.15)

axs = axs.flatten()
for p in range(4):
	#plt.figure()
	#ax=plt.subplot(p)
	file ='porosity_table_' + files[p]
	data =   np.loadtxt(open(file  + '.csv' , 'rb'), delimiter=',')
	#t = np.linspace(0,data.shape[0], data.shape[0])
	
	markers = ["o", "v", "^", "<", ">", "8"]	
	colors = plt.cm.Greys(np.linspace(0.5,1,len(t)))
	#x = list(reader)
	#result = numpy.array(x).astype("float")
	for i in range(len(t)):

		#label = "timestep " + str(t[i] +1)
		#plt.plot(dist, data[t[i]], label = labels[i] )
		axs[p].plot(dist, data[t[i]], label = labels[i], color = colors[i], marker=markers[i], linestyle='solid' )# color='green',linewidth=2, markersize=12)


	axs[p].set_ylim(0,1.1)
	axs[p].set_xlim(0, 150)

	axs[p].set_xlabel('Distance from root surface [$\mu$m] ')
	axs[p].set_ylabel('Porosity [%]')
	axs[p].title.set_text(names[p])
	handles, labels = axs[p].get_legend_handles_labels()
fig.legend(handles, labels,  ncol=6,loc='lower center')#loc='upper center',ncol=6,labelspacing=0.,borderaxespad=0.)
	# plt.legend()
fig.savefig('Fig.5' + '.png', bbox_inches='tight')
fig.show()