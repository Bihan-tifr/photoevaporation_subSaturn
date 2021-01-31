#This is a python program for photoevaporation
#This code can be used in other codes
#If initial mass radius period samples are given, it can compute final radius and mass samples

import numpy as np
import matplotlib.pyplot as plt
import warnings
from astropy import constants as const
import warnings
from M_R import ck17_M_to_R as Rck
from M_R import gen_R as Rgn
from M_R import lopez_R as Rlf

warnings.filterwarnings("ignore")
#constants
pi=np.pi
M_earth=1	#unit mass
R_earth=1	#unit length	
Myr=1	#unit time
M_sun=3*10**5	#in earth mass	
G=1.51*10**21	#in proper units
L0=4.86*10**(28)	#in proper units
au=2.34*10**4  #earth radius
n=0.1


#cgs values(used to calculate beta)
me_cgs=const.M_earth.cgs.value
re_cgs=const.R_earth.cgs.value
G_cgs=const.G.cgs.value
myr_cgs=10**6*(np.pi*10**7)
def P_a(Ms,P):
	P=P*24*3600
	Ms0=Ms*me_cgs
	a=(G_cgs*Ms0*(P/2*pi)**2)**(1/3)
	return a/re_cgs		#distance in earth radius unit
def lumX(t,Ms):
	t_sat=300
	L_sat=10**(-3.5)*L0
	if(t<=t_sat):
		return L_sat
	if(t>t_sat):
		return L_sat*(t/t_sat)**(-1.5)
def XUV(Lx):
	Lxc=Lx*((10**33)/L0)
	log_L_euv = 0.860 * np.log10(Lxc) +4.8
	Leuv = 10.**log_L_euv
	log_L_xuv = np.log10(Lxc + Leuv)
	return 10**(log_L_xuv)*(L0*10**(-33))
	
def kfn(Mp,Ms,Rp,a):
	R_roche=a*(Mp/3*Ms)**(1/3)
	return 1-(3/2)*(Rp/R_roche)+(1/2)*(Rp/R_roche)**3
		
def Mdot(t,Mp,Rp,a,Ms):
	Lx=lumX(t,Ms)
	Lxuv=XUV(Lx)
	Fxuv=Lxuv/(4*np.pi*a**2)
	rho=Mp/((4/3)*np.pi*Rp**3)
	b=1.5
	k=kfn(Mp,Ms,Rp,a)
	Mp1=n*(3*b**2*Fxuv)/(4*G*k*rho)			
	return Mp1
	
def photoevolve(Msamples,Psamples,Rsamples,M_s,t_final,dt,par=None):
	if(par==None):
		Ms=M_s*M_sun
		alist=P_a(Ms,Psamples)
		
		
		M=Msamples
		a=alist
		p=Psamples
		R=Rsamples
		t=0
		while(t<t_final):
			k1=dt*Mdot(t,M,Rck(M),a,Ms)
			k2=dt*Mdot(t+dt/2,M-k1,Rck(M-k1),a,Ms)
			k3=dt*Mdot(t+dt/2,M-k2,Rck(M-k2),a,Ms)
			k4=dt*Mdot(t+dt,M-k3,Rck(M-k3),a,Ms)
			
			M=M-k1/6-k2/3-k3/3-k4/6
			R=Rck(M)
			t+=dt
			
		R_final=np.array(R)
			
		return np.array(R_final)
	else:
		Ms=M_s*M_sun
		alist=P_a(Ms,Psamples)
		
		
		M=Msamples
		a=alist
		p=Psamples
		R=Rsamples
		t=0
		while(t<t_final):
			k1=dt*Mdot(t,M,Rgn(t,M,par),a,Ms)
			k2=dt*Mdot(t+dt/2,M-k1,Rgn(t+dt/2,M-k1,par),a,Ms)
			k3=dt*Mdot(t+dt/2,M-k2,Rgn(t+dt/2,M-k2,par),a,Ms)
			k4=dt*Mdot(t+dt,M-k3,Rgn(t+dt,M-k3,par),a,Ms)
			
			M=M-k1/6-k2/3-k3/3-k4/6
			R=Rgn(t+dt,M,par)
			t+=dt
			
		R_final=np.array(R)
		M_final=np.array(M)
			
		return [R_final,M_final]		
		
def photoevolve_L(Msamples,M_core,Psamples,M_s,t_initial,t_final,dt):
	Ms=M_s*M_sun
	a=P_a(Ms,Psamples)
	
	Mp=Msamples
	t=t_initial
	while(t<t_final):
		k1=dt*Mdot(t,Mp,Rlf(t,Mp,M_core,a),a,Ms)
		k2=dt*Mdot(t+dt/2,Mp-k1,Rlf(t+dt/2,Mp-k1,M_core,a),a,Ms)
		k3=dt*Mdot(t+dt/2,Mp-k2,Rlf(t+dt/2,Mp-k2,M_core,a),a,Ms)
		k4=dt*Mdot(t+dt,Mp-k3,Rlf(t+dt,Mp-k3,M_core,a),a,Ms)
		
		Mp=Mp-k1/6-k2/3-k3/3-k4/6
		R=Rlf(t+dt,Mp,M_core,a)
		t+=dt
		
	R_final=np.array(R)
	M_final=np.array(Mp)
		
	return [R_final,M_final]
		
		
		
		
		









