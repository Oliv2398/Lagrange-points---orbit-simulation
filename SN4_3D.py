##### Olivier Roth 
##### SN4 
##### Equation trajectoire M (en 3D)



### Importation des modules

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

tmp_1=time.time()



### Constantes

G=6.67259e-11 # SI
M_S=1.989e30 # kg
M_T=5.972e24 # kg
dS_T=149.597870700e9 # m
d1=1.4917449889577318e9 # m, distance Terre - point L1
d2=-1.501728936924601e9 # m, distance Terre - point L2
m= 1800  # kg, masse de l'objet
v0=np.sqrt(G*M_S/dS_T) # = 29785.294991204377 m/s , vitesse pour la Terre
v0m1=np.sqrt(G*(M_S+M_T)/(dS_T+d1)) #  vitesse au point L1
v0m2=np.sqrt(G*(M_S+M_T)/(dS_T+d2)) # vitesse au point L2
v1=-149.175379 # correction de vitesse pour L1
v2=147.84918 # correction de vitesse pour L2



### Fonctions

fS_Tx= lambda x_t, y_t, z_t : -G*M_T*M_S*x_t/((x_t**2+y_t**2+z_t**2)**(3/2)) # Force du Soleil sur la Terre en x
fS_Ty= lambda x_t, y_t, z_t : -G*M_T*M_S*y_t/((x_t**2+y_t**2+z_t**2)**(3/2)) # Force du Soleil sur la Terre en y
fS_Tz= lambda x_t, y_t, z_t : -G*M_T*M_S*z_t/((x_t**2+y_t**2+z_t**2)**(3/2)) # Force du Soleil sur la Terre en z


fx = lambda x, y, z, X_T, Y_T, Z_T : -G*m*M_S*x/((x**2+y**2+z**2)**(3/2)) - G*m*M_T*(x-X_T)/(((x-X_T)**2+(y-Y_T)**2+(z-Z_T)**2)**(3/2)) # force en x exercée sur un objet de masse m
fy = lambda x, y, z, X_T, Y_T, Z_T : -G*m*M_S*y/((x**2+y**2+z**2)**(3/2)) - G*m*M_T*(y-Y_T)/(((x-X_T)**2+(y-Y_T)**2+(z-Z_T)**2)**(3/2)) # force en y exercée sur un objet de masse m
fz = lambda x, y, z, X_T, Y_T, Z_T : -G*m*M_S*z/((x**2+y**2+z**2)**(3/2)) - G*m*M_T*(z-Z_T)/(((x-X_T)**2+(y-Y_T)**2+(z-Z_T)**2)**(3/2)) # force en z exercée sur un objet de masse m



def Euler_Terre(x0, y0, v0) :
    ttab=np.zeros(N)
    xtab=np.zeros(N)
    ytab=np.zeros(N)
    ztab=np.zeros(N)
    vxtab=np.zeros(N)
    vytab=np.zeros(N)
    vztab=np.zeros(N)
    xtab[0]=x0
    ytab[0]=y0
    ztab[0]=0
    vxtab[0]=0
    vytab[0]=v0
    vztab[0]=0
    ttab[0]=0
    for i in range(N-1) :
    	ttab[i+1] = ttab[i]+dt
    	xtab[i+1] = xtab[i]+dt*vxtab[i]
    	ytab[i+1] = ytab[i]+dt*vytab[i]
    	ztab[i+1] = ztab[i]+dt*vztab[i]
    	vxtab[i+1] = vxtab[i]+dt*fS_Tx(xtab[i],ytab[i],ztab[i])/M_T
    	vytab[i+1] = vytab[i]+dt*fS_Ty(xtab[i],ytab[i],ztab[i])/M_T
    	vztab[i+1] = vztab[i]+dt*fS_Tz(xtab[i],ytab[i],ztab[i])/M_T
    return xtab, ytab, ztab


def Euler_objet(d0, d1, v0, y0, z0, vz0) :
    ttab=np.zeros(N)
    xtab=np.zeros(N)
    ytab=np.zeros(N)
    ztab=np.zeros(N)
    vxtab=np.zeros(N)
    vytab=np.zeros(N)
    vztab=np.zeros(N)
    xtab[0]=d0-d1
    ytab[0]=y0
    ztab[0]=z0
    vxtab[0]=0
    vytab[0]=v0
    vztab[0]=vz0
    ttab[0]=0
    for i in range(N-1) :
    	ttab[i+1] = ttab[i]+dt
    	xtab[i+1] = xtab[i]+dt*vxtab[i]
    	ytab[i+1] = ytab[i]+dt*vytab[i]
    	ztab[i+1] = ztab[i]+dt*vztab[i]    	
    	vxtab[i+1] = vxtab[i]+dt*(fx(xtab[i],ytab[i],ztab[i], X_T[i], Y_T[i],Z_T[i]))/m
    	vytab[i+1] = vytab[i]+dt*(fy(xtab[i],ytab[i],ztab[i], X_T[i], Y_T[i],Z_T[i]))/m
    	vztab[i+1] = vztab[i]+dt*(fz(xtab[i],ytab[i],ztab[i], X_T[i], Y_T[i],Z_T[i]))/m
    return xtab, ytab, ztab



def plot_3D(x, y, z, xt, yt, zt, b) : # pour un graph en 3D
    ax = plt.axes(projection='3d')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.plot(x, y, z, '--r', label="object in L"+b)
    ax.plot(xt, yt, zt, '-b', label="Earth")		
    #plt.title("Trajectoire de la Terre et d'un objet au point L"+b)
    plt.legend()
    plt.show()


def plot_anim(N, x, y) : # pour un plot animé
    plt.ion()
    n=0
    while n<N:
	    plt.axis([-2e11,2e11,-5e4,5e4])
	    plt.plot(x[n],y[n],'xr')
	    plt.draw()
	    plt.pause(1e-2)
	    n+=5000
    plt.ioff()
    plt.show()	




### Programme principal

dt=100 # s
N=315582 # *dt -> nombre de secondes dans une année (une révolution de la Terre), 365.256363051 est le nombre de jours dans une année.



X_T, Y_T, Z_T=Euler_Terre(dS_T, 0, v0) # pour trouver les positions de la Terre


#X_M, Y_M, Z_M=Euler_objet(dS_T, d1, v0m1+v1, 0, 0, 0) # objet en L1
#X_M2, Y_M2, Z_M2=Euler_objet(dS_T, d2, v0m2+v2, 0, 0, 0) # objet en L2

X_M3, Y_M3, Z_M3=Euler_objet(dS_T, d1, v0m1+v1, 0, 1e4, -0.01) #objet en orbite autour de L1





tmp_2=time.time()

tmp=tmp_2-tmp_1
print(tmp, "secondes pour exécuter le programme")






### Plot


#plot_3D(X_M, Y_M, Z_M, X_T, Y_T, Z_T, "1") # graph 3D pour L1
#plot_3D(X_M2, Y_M2, Z_M2, X_T, Y_T, Z_T, "2") # graph 3D pour L2


"""
plt.plot(X_M3, Z_M3, '-r', label='trajectoire de type Lissajous')
plt.plot(X_T, Z_T,'-b', label='trajectoire de la Terre ')
plt.title("Coupe selon le plan x0z d'une orbite de Lissajous")
plt.xlabel('x')
plt.ylabel('z')
plt.legend()
plt.show()
"""

#plot_anim(N, Y_M3, Z_M3)


plot_3D(X_M3, Y_M3, Z_M3, X_T, Y_T, Z_T, "1") # orbite autour de L1


"""
X_M4, Y_M4, Z_M4=Euler_objet(dS_T, d1, v0m1+v1, dt, N, 0, 10, -1.8) 
#plot_3D(X_M4, Y_M4, Z_M4, X_T, Y_T, Z_T, "1")
plt.plot(X_M4, Z_M4)
plt.show()
"""


