import numpy as np
import rebound
import bin_delete


N=100
N2=int(1.5*N)

##N is the true number of particles to delete
##add particles as usual...and initialize simulation as usual.
#sim=...
bin_delete.delete_bins(sim)

to_del=range(N+1, len(sim.particles))[::-1]
for idx in to_del:
	# print(type(idx))
	sim.remove(index=idx)
print(len(sim.particles))