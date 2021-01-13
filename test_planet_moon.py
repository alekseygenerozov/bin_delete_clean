from bash_command import bash_command as bc

import sys
# sys.path.append('/usr/local/lib/python2.7/dist-packages/')
import rebound
import bin_delete
import numpy as np



a_bin=3.0e-3
inc_bin=np.pi/3.
e_bin=0.7
tol=1.0e-5

sim2 = rebound.Simulation()
sim2.G = 1.
sim2.add(m = 1.) # Star
sim2.add(m = 3e-6, a=1., inc=0.) #Planet
sim2.add(m = 3e-10, a=a_bin, e=e_bin, inc=inc_bin, primary=sim2.particles[1]) #Moon

sim3 = rebound.Simulation()
sim3.G = 1.
sim3.add(m = 1.) # Star
for i in range(1, 4):
	sim3.add(m = 3e-6, a=i*2, primary=sim3.particles[0]) #Planets
sim3.add(m=3e-10, a=a_bin, e=e_bin, inc=inc_bin, primary=sim3.particles[3])

def test_bin_props():
	x=bin_delete.bin_props(sim2.particles[1], sim2.particles[2])
	assert abs(x[1]-a_bin)/a_bin<tol
	assert abs(x[2]**0.5-e_bin)/e_bin<tol
	print(x[-2]*np.pi/180.)
	assert abs(x[-2]*np.pi/180.-inc_bin)/inc_bin<tol

def test_bin_find3():
	bins=bin_delete.bin_find_sim(sim3)
	x=bin_delete.bin_props(sim3.particles[3], sim3.particles[4])
	assert abs(x[1]-a_bin)/a_bin<tol
	assert (len(bins)==1)
	assert set(bins[0,[1,2]])=={3,4}

def test_bin_delete():
	bin_delete.delete_bins(sim3)
	assert len(sim3.particles)==4
	masses=np.array([pp.m for pp in sim3.particles[1:]])
	assert(np.all(masses==np.array([3e-6, 3e-6, 3e-10])))


