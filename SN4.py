##### Olivier Roth 
##### SN4 
##### Positions L1, L2 par recherche de zéros 


### Importation des modules

import numpy as np
import matplotlib.pyplot as plt
import math
import time

tmp_1=time.time()


### Fonctions

def f1(x):
	return -1/((1-x)**2) + mu/(x**2) +1 -(1+mu)*x

def fp1(x):
	return -2*x/((1-x)**4) -2*mu/(x**3) -1+mu

def f2(x):
	return -1/((1-x)**2) - mu/(x**2) +1 -(1+mu)*x

def fp2(x):
	return -2*x/((1-x)**4) + 2*mu/(x**3) -1+mu


def dicho(f, a, b, eps):
	n=0
	while (b-a) > eps :
		c= (a+b)/2
		if f(a)*f(c)<=0 :
			b=c
		else :
			a=c
	return c

def newr(f,fp,n,x0):
	i=0
	u=x0
	while i<n:
		u=u-f(u)/fp(u)
		i+=1
	return u

### Programme principal

mu=1/332830
eps=0.0000001


print("Positions des points L1 et L2 par recherche de zéros :\n\n")

print("Distance Terre-pt L1 (par dichotomie) :", dicho(f1,0.001,1,eps)*149.6e6, "km", "soit", np.round(dicho(f1,0.001,1,eps),5),"ua de la Terre." )
print("\nDistance Terre-pt L1 (par Newton-Raphson) :",newr(f1,fp1,100,0.001)*149.6e6, "km", "soit", np.round(newr(f1,fp1,100,0.001),5),"ua de la Terre.\n")


print("\nDistance Terre-pt L2 (par dichotomie) :", dicho(f2,-1+eps,0-eps, eps)*149.6e6, "km", "soit", np.round(dicho(f2,-1+eps,0-eps, eps),5), "ua de la Terre.")
print("\nDistance Terre-pt L2 (par Newton-Raphson) :",newr(f2,fp2,100,0.001)*149.6e6, "km", "soit", np.round(newr(f2,fp2,100,0.001),5), "ua de la Terre.\n\n")


print("~~~~~~~~~~~~~~")

tmp_2=time.time()

tmp=tmp_2-tmp_1
print('\n', np.round(tmp,5), "secondes pour exécuter le programme.")

