#Aim is to constrain a,b values of toy model to successfully create sub saturn desert.
#Draw period samples from the distribution/observed value
#for each period value randomly select a logM value from uniform distribution.(1<M<200): so (0<logM<2.3)
#find radius samples using "R=cM^a(t)" relationship
#find proper "a" value which successfully creates the desert. To do this, select a range of parameters and separately work out for each and pick the best. 
#how to define the best? (1) number of planets in the desert is very small (2)minimum xi square

run=7
import numpy as np
import matplotlib.pyplot as plt
import warnings
from astropy import constants as const
from scipy.optimize import minimize
import warnings
from scipy.optimize import fsolve


#Mazeh boundaries:

log_p=np.linspace(-0.1,0.9,1000)
logR_up =-0.33*log_p+1.17
logR_dn = 0.68*log_p+0.3

#print(M)

import sys
sys.path.append('/home/bihan/Documents/project codes/dp')
from M_R import gen_R as Rgn
from photoevolve import photoevolve


import pandas as pd
df=pd.read_csv('Planets.csv')
porbd=df['pl_orbper']
porb=porbd.to_numpy()
orbp_s=porb[(porb>0.7)*(porb<8)]
p_s=np.random.choice(orbp_s,2000)

t_final=5000  # in Myr
dt= 0.1	#in Myr
M_s=1		#in solar mass

Mlist=[]
plist=[]
Rlist=[]
par=[0.9,0.9]
for i in range(2000):
	p=p_s[i]
	logm=np.random.uniform(0,0.68*np.log10(p)+1.2)
	M=10**(logm)
	R=Rgn(0,M,par)
	ub=-0.33*np.log10(p)+1.17

	plist.append(p)
	Mlist.append(M)
	Rlist.append(R)

Psamples=np.array(plist)
Msamples=np.array(Mlist)
Rsamples=np.array(Rlist)

R_final,M_final=photoevolve(Msamples,Psamples,Rsamples,M_s,t_final,dt,par)

#initial R-p
plt.figure(1)
plt.xlabel("Period in days")
plt.ylabel("Radius in earth radius")
plt.plot(10**(log_p),10**(logR_up),'k',label="desert boundary")
plt.plot(10**(log_p),10**(logR_dn),'k')
plt.scatter(Psamples,Rsamples,s=2,label="Newborn planets: param={}".format(par))
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend(loc='best')
plt.savefig("gm_init_RP{}.png".format(run))

#final R-P
plt.figure(2)
plt.xlabel("Period in days")
plt.ylabel("Radius in earth radius")
plt.plot(10**(log_p),10**(logR_up),'k',label="desert boundary")
plt.plot(10**(log_p),10**(logR_dn),'k')
plt.scatter(Psamples,R_final,s=2,label="Evolved planets: param={}".format(par))
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend(loc='best')
plt.savefig("gm_final_RP{}.png".format(run))

#Initial M-P relation
plt.figure(3)
plt.xlabel("Period in days")
plt.ylabel("Mass in earth mass")
plt.scatter(Psamples,Msamples,s=2,label="Newborn planets: param={}".format(par))
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend(loc='best')
plt.savefig("gm_init_MP{}.png".format(run))

#Final M-P relation
plt.figure(4)
plt.xlabel("Period in days")
plt.ylabel("Mass in earth mass")
plt.scatter(Psamples,M_final,s=2,label="Evolved planets: param={}".format(par))
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend(loc='best')
plt.savefig("gm_final_MP{}.png".format(run))

#Initial M-R relation
plt.figure(5)
plt.xlabel(r'Mass in $M_\oplus$',fontsize=11)
plt.ylabel(r'radius in $R_\oplus$',fontsize=11)
plt.scatter(Msamples,Rsamples,s=2,label="Newborn planets, param={}".format(par))
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend(loc='best')
plt.savefig("gm_initial_MR{}.png".format(run))

#Final M-R relation
plt.figure(6)
plt.xlabel(r'Mass in $M_\oplus$',fontsize=11)
plt.ylabel(r'radius in $R_\oplus$',fontsize=11)
plt.scatter(M_final,R_final,s=2,label="Evolved planets, param={}".format(par))
plt.xscale('log')
plt.yscale('log')
plt.grid(True)
plt.legend(loc='best')
plt.savefig("gm_final_MR{}.png".format(run))



'''
#period samples
import pandas as pd
df=pd.read_csv('Planets.csv')
porbd=df['pl_orbper']
porb=porbd.to_numpy()
Mlist=[]
plist=[]
Rlist=[]
par=[0.6,0.8]
import sys
sys.path.append('/home/bihan/Documents/project codes/dp')
from M_R import gen_R as Rgn

for i in range(len(porb)):
	logp=np.log10(porb[i])

	u_bound=-0.33*logp+1.17
	l_bound=0.68*logp
	Ru=10**(u_bound)
	Rl=10**(l_bound)
	if(u_bound>l_bound):
		logm=np.random.uniform(-0.5,2.3)
		M=10**(logm)
		R_init=Rgn(0,M,par)
		if(R_init<=5):
			Mlist.append(M)
			Rlist.append(R_init)
			p=porb[i]
			plist.append(p)


p_samples=np.array(plist)
M_samples=np.array(Mlist)
R_samples=np.array(Rlist)

plt.scatter(p_samples,R_samples,s=2)
plt.figure(2)
plt.scatter(p_samples,M_samples,s=2)
plt.xlabel("period")
plt.ylabel("initial mass")
plt.xscale('log')
plt.yscale('log')
plt.grid('True')
plt.figure(3)
plt.scatter(M_samples,R_samples,s=2)
plt.xlabel("initial mass")
plt.ylabel("initial radius")
plt.xscale('log')
plt.yscale('log')
plt.grid('True')
plt.show()

t_final=5000  # in Myr
dt= 0.1	#in Myr
M_s=1		#in solar mass

from photoevolve import photoevolve

R_final,M_final= photoevolve(M_samples,p_samples,R_samples,M_s,t_final,dt,par)


plt.figure(4)
plt.scatter(p_samples,R_final,s=2,label="evolved planets")
plt.plot(10**(log_p),10**(logR_up),'k',label="desert boundary")
plt.plot(10**(log_p),10**(logR_dn),'k')
plt.xscale('log')
plt.yscale('log')
plt.grid('True')
plt.legend(loc='best')
plt.xlabel("Period in days")
plt.ylabel("Final radius")
plt.figure(5)
plt.scatter(p_samples,M_final,s=2)
plt.xlabel("period")
plt.ylabel("Final mass")
plt.xscale('log')
plt.yscale('log')
plt.grid('True')
plt.figure(6)
plt.scatter(M_final,R_final,s=2)
plt.xlabel("Final mass")
plt.ylabel("Final radius")
plt.xscale('log')
plt.yscale('log')
plt.grid('True')


plt.show()
'''
