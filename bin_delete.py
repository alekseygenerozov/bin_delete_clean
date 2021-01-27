import rebound
import numpy as np
from itertools import combinations


def bin_find_sim(sim):
	##Ensure we are in the com frame of the simulation.
	sim.move_to_com()
	##Particles in simulation
	ps = sim.particles
	##mass of of primary 
	m0 = ps[0].m
	bin_indics=[]
	for i1, i2 in combinations(list(range(1, sim.N)),2): # get all pairs of indices/start at 1 to exclude SMBH
		oo = ps[i1].calculate_orbit(primary=ps[i2])
		a_bin=oo.a
		e_bin=oo.e
		m1,m2 =ps[i1].m, ps[i2].m
		com =(m1*ps[i1]+m2*ps[i2])/(m1+m2)
		com_d = (com.x**2.+com.y**2.+com.z**2.)**0.5
		##*approximate* Hill sphere condition. 
		inside_hill=(a_bin<((m1+m2)/m0)**(1./3.)*com_d)
		dp = ps[i1] - ps[i2]
		d2 = dp.x*dp.x+dp.y*dp.y+dp.z*dp.z

		##If the kinetic energy is less than the potential energy 
		if ((a_bin>0) and (inside_hill)):
			rh=(((m1+m2)/m0)**(1./3.)*com_d)
			vh=rh*(m0/com_d**3.)**0.5
			bin_indics.append([sim.t, i1, i2, d2**0.5, a_bin, a_bin/(((m1+m2)/m0)**(1./3.)*com_d), e_bin, rh, vh])

	return np.array(bin_indics)

def delete_bins(sim):
	'''
	Delete all binaries from a simulation

	Input--rebound simulation "sim"
	'''	

	sim.move_to_com()
	##Look for binaries
	bins=bin_find_sim(sim)
	##Loop until all binaries are deleted. Doing a loop is only really necessary with tidal criterion which is not really used right now
	while len(bins)>0:
		##Delete in reverse order (else the indices would become messed up)
		to_del=(np.sort(np.unique(bins[:,1]))[::-1]).astype(int)
		##Delete all binaries.
		for idx in to_del:
			sim.remove(index=int(idx))
		bins=bin_find_sim(sim)
