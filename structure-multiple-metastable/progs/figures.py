
import matplotlib

import matplotlib.pyplot as plt
#import matplotlib.pylab 
import matplotlib.mlab as mlab
import numpy as np

# USE TEX:
#matplotlib.use('ps')
matplotlib.rc('text', usetex=True)



	


base = "ising_4x4_beta_1.0"


fig = plt.figure(figsize=(7, 2))
#ax = fig.add_subplot(111)



#plt.xticks(np.arange(0, 70000, 10000))

#plt.yticks(np.linspace(0, 1, 5))



threshold = 1e-2


plt.axes([0.3,0.2,0.95-0.3,0.95-0.2])

files = [ 
			["ising_4x4_beta_1.0", 1], 
			  #["potts_fast_slow_3x4_beta_1.2", 2],
			  ["potts_fast_slow_3x3_beta_1.2", 2],
			  #["potts_slow_fast_3x4_beta_1.2", 2],
			  ["potts_slow_fast_3x3_beta_1.2", 2],
			  #["potts_fast_slow_3x4_beta_1.2", 2]
			  ["potts_succession_q4_3x3_beta_1.5", 3]
			]

for filename, num in files:
	
	
	for kind in ['P', 'C']:
		
		
		plt.clf()
		plt.axes([0.3,0.2,0.95-0.3,0.95-0.2])
		
		for i in range(num+1):
			
			
			
			if kind=='C' and i==0:
				continue
			
			if kind=='C' and len(data)>66000:
				continue
			
			if kind=='C':
				plt.subplot(1, num, i)
			else:
				plt.subplot(1, num+1, i+1)
			
			data = mlab.load("%s_%s%d.dat" % (filename, kind, i))
			
			num_elems = len(data)
			
			t = np.arange(num_elems)
			
			
			
			plt.hold(True)
			
			max_val = max(abs(data)) * threshold
			
			if kind == 'P':
				size=4
				width = 0.2
			else:
				size=1
				width = 0.1
	
	
			
			
			
			
			#plt.grid(linestyle=':', color='0.3', linewidth=0.1, alpha=0.3)
			
			
			plt.plot(t[abs(data)>max_val], data[abs(data)>max_val], 'o', markersize=size, markeredgewidth=width,
						scalex = False, scaley=True, markeredgecolor='b')
			#plt.axis((-1,66000.,-0.01,1.01))
			plt.xlim((-0.05*num_elems,1.05*num_elems))
				
				
				
			tick_pos = np.linspace(0,1,num+1)*num_elems
			
			# vertical grid lines
			for j in tick_pos:
				plt.axvline(x=j, linestyle=':', color='0.3', linewidth=0.1, alpha=0.3)
				
			plt.axhline(0, linestyle='-', color='0.3', linewidth=0.1, alpha=0.3)
			
			
			plt.xlabel('$\sigma$')
			plt.ylabel('$%s_%d(\sigma)$' % (kind, i) )
			
			
			
			
			labels = []
			for j in range(num+1):
				labels.append( '$\mathbf{%d}$' % j )
				
				
				plt.xticks(tick_pos, labels)
				

			
		outfile_name = '%s_%s.eps' % (filename, kind)
		print "Creating ", outfile_name
		plt.savefig(outfile_name)
			
			#plt.clf()
			
			#plt.axes([0.3,0.2,0.95-0.3,0.95-0.2])
			
	#ax.plot(x, y)

##plt.show()


