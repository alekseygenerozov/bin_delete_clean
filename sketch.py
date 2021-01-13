import numpy as np
import rebound
import bin_delete


##N is the true number of stars you want in your simulation
N=100
##Add a buffer
N2=int(1.5*N)

##add particles as usual...and initialize simulation as usual.
#sim=...
##Delete binaries
bin_delete.delete_bins(sim)

#Delete excess particles so we are left with N stars + (central object.)
to_del=range(N+1, len(sim.particles))[::-1]
for idx in to_del:
	# print(type(idx))
	sim.remove(index=idx)
print(len(sim.particles))