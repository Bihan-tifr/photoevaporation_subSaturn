#this code explores various mass radius relationships. 
#(1) Chen Kipping (2) Lopez-Fortney (3) Generic Models

import numpy as np

def ck17_M_to_R(M):
	if(type(M)==np.float64):
		if(M<=2):
			R=1.01*M**(0.28)
		else:
			R=0.82*M**(0.59)
		return R
		
	elif(type(M)==np.ndarray):
		Rlist=[]
		for i in range(len(M)):
			m=M[i]
			if(m<=2):
				R=1.01*m**(0.28)
			else:
				R=0.82*m**(0.59)	
			Rlist.append(R)
			
		return np.array(Rlist)

		
def lopez_R(t,Mp,M_core,a):
	au=2.34*10**4  #earth radius
	R_earth=1
	M_earth=1
	if(type(Mp)==np.ndarray):
		I=np.ones(len(Mp))
		fenv=(I-M_core/Mp)
		R_core=(M_core)**(0.25)
		flr=0.934*(au/(a))**2
		Renv= 2.06*R_earth*(Mp/M_earth)**(-0.21)*(fenv/0.05)**(0.59)*flr**(0.044)*(t/5000)**(-0.11)
		Rp= Renv+R_core
		return Rp
	else:
		fenv=(1-M_core/Mp)
		flr=0.934*(au/(a))**2
		Renv= 2.06*R_earth*(Mp/M_earth)**(-0.21)*(fenv/0.05)**(0.59)*flr**(0.044)*(t/5000)**(-0.11)
		Rp= Renv+R_core
		return Rp

def gen_R(t,M,par):
	a,b=par
	if(type(M)==np.ndarray):
		Rlist=[]
		for i in range(len(M)):
			m=M[i]
			if(m<=2):
				R=1.01*m**(a+(0.28-a)*(t/5000))
			else:
				R=0.82*m**(b+(0.59-b)*(t/5000))	
			Rlist.append(R)
			
		return np.array(Rlist)	
	else:
		if(M<=2):
			R=1.01*M**(a+(0.28-a)*(t/5000))
			return R
		else:
			R=0.82*M**(b+(0.59-b)*(t/5000))
			return R	

