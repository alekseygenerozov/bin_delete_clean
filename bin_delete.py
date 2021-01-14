import rebound
import numpy as np
from itertools import combinations


def bin_props(p1, p2):
	'''
	Auxiliary function to get binary properties for two particles. 

	p1 and p2 -- Two particles from a rebound simulation.
	'''
	dp = p1 - p2   
	d2 = dp.x*dp.x+dp.y*dp.y+dp.z*dp.z
	# print(d2)

	##Calculate com velocity of two particles...
	##Masses
	m1=p1.m
	m2=p2.m
	com = (m1*p1+m2*p2)/(m1+m2)
	##Particle velocities in com frame
	p1_com = p1 - com
	p2_com = p2 - com
	v12 = (p1_com.vx**2.+p1_com.vy**2.+p1_com.vz**2.)
	v22 = (p2_com.vx**2.+p2_com.vy**2.+p2_com.vz**2.)
	##Difference in the forces acting on the two particles;
	ft = np.array([m2*(p2.ax)-m2*(com.ax), m2*(p2.ay)-m2*(com.ay), m2*(p2.az)-m2*com.az])
	##Unit vector pointing from particle 2 to particle 1
	rhat = np.array(dp.xyz)/d2**0.5
	f12 = m1*m2/d2*rhat 
	##Tidal force that star 2 experiences (NOT CURRENTLY USED FOR BINARY DELETION)
	ft = ft - f12
	ft = np.linalg.norm(ft)

	##Kinetic and potential energies
	ke = 0.5*m1*v12+0.5*m2*v22
	##Potential energy; Assumes G = 1
	pe = (m1*m2)/d2**0.5

	##Distance of binary center of mass from COM of system (should be near central SMBH)
	com_d=(com.x**2.+com.y**2.+com.z**2.)**0.5
	a_bin=(m1*m2)/(2.*(pe-ke))
	##Angular momentum in binary com
	j_bin=m1*np.cross(p1_com.xyz, p1_com.vxyz)+m2*np.cross(p2_com.xyz, p2_com.vxyz)
	##Angular momentum of binary com
	j_com=(m1+m2)*np.cross(com.xyz, com.vxyz)

	#Inclination of star's orbit wrt the binary's orbit win the disk *in degrees*
	inc=np.arccos(np.dot(j_bin, j_com)/np.linalg.norm(j_bin)/np.linalg.norm(j_com))*180./np.pi
	#inc=np.nan
	mu=m1*m2/(m1+m2)
	##Eccentricity of the binary *squared*
	e_bin=(1.-np.linalg.norm(j_bin)**2./((m1+m2)*a_bin)/(mu**2.))

	return com_d, a_bin, e_bin, p1_com, p2_com, d2, inc, ft

 
def bin_find_sim(sim):
	'''
	Find indices of all binaries in a rebound simulation.


	Input--rebound simulation "sim"

	'''
	##Ensure we are in the com frame of the simulation.
	sim.move_to_com()
	##Integrate forward for a small time to ensure that the accelerations are initialized (not necessary since we are not using tidal condition.)
	# sim.integrate(sim.t+sim.t*1.0e-14)
	
	ps = sim.particles
	##mass of of primary 
	m0 = ps[0].m
	bin_indics=[]
	for i1, i2 in combinations(list(range(1, sim.N)),2): # get all pairs of indices/start at 1 to exclude SMBH
		com_d, a_bin, e_bin, p1_com, p2_com, d2, inc, ft = bin_props(ps[i1], ps[i2])
		m1,m2 =ps[i1].m, ps[i2].m
		##Hill sphere condition.
		inside_hill=(a_bin<((m1+m2)/m0)**(1./3.)*com_d)
		##Tidal condition not currently used...
		# tidal_2 = (m1*m2/d2>ft)

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
