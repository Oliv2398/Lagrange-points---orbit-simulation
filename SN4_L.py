##### Olivier Roth 
##### SN4 
##### Trajectoire Lune


### Importation des modules

import numpy as np
import matplotlib.pyplot as plt
import time


### Fonctions

tmp_1=time.time()


fS_Tx= lambda x_t, y_t : -G*M_T*M_S*x_t/((x_t**2+y_t**2)**(3/2)) # Force du Soleil sur la Terre en x
fS_Ty= lambda x_t, y_t : -G*M_T*M_S*y_t/((x_t**2+y_t**2)**(3/2)) # Force du Soleil sur la Terre en y

fT_Lx= lambda x, y, x_t, y_t : -G*M_T*M_L*(x-x_t)/(((x-x_t)**2+(y-y_t)**2)**(3/2))-G*M_L*M_S*x/((x**2+y**2)**(3/2)) # Force sur la Lune en x
fT_Ly= lambda x, y, x_t, y_t : -G*M_T*M_L*(y-y_t)/(((x-x_t)**2+(y-y_t)**2)**(3/2))-G*M_L*M_S*y/((x**2+y**2)**(3/2)) # Force sur la Lune en y



def Euler(fT_Lx, fT_Ly, fS_Tx ,fS_Ty, dS_T, d0, v0, v1, dt, N) :
    ttab=np.zeros(N)
    xtab=np.zeros(N)
    x2tab=np.zeros(N)
    ytab=np.zeros(N)
    y2tab=np.zeros(N)
    vxtab=np.zeros(N)
    vx2tab=np.zeros(N)
    vytab=np.zeros(N)
    vy2tab=np.zeros(N)
    
    xtab[0]=dS_T
    x2tab[0]=d0
    vytab[0]=v0
    vy2tab[0]=v1
    
    for i in range(N-1) :
    	ttab[i+1] = ttab[i]+dt
    	xtab[i+1] = xtab[i]+dt*vxtab[i]
    	x2tab[i+1] = x2tab[i]+dt*vx2tab[i]
    	ytab[i+1] = ytab[i]+dt*vytab[i]
    	y2tab[i+1] = y2tab[i]+dt*vy2tab[i]
    	vxtab[i+1] = vxtab[i]+dt*fS_Tx(xtab[i],ytab[i])/M_T
    	vx2tab[i+1]=vx2tab[i]+dt*fT_Lx(x2tab[i],y2tab[i], xtab[i], ytab[i])/M_L
    	vytab[i+1] = vytab[i]+dt*fS_Ty(xtab[i],ytab[i])/M_T
    	vy2tab[i+1] = vy2tab[i]+dt*fT_Ly(x2tab[i],y2tab[i], xtab[i], ytab[i])/M_L
    return ttab, xtab, ytab, x2tab, y2tab


def plot_traj(x, y, xt, yt, a, b) : # plot des résultats des équations différentielles
    plt.axis([-a*1e11,a*1e11,-a*1e11,a*1e11])
    plt.axis('equal')
    plt.plot(xt, yt, '-b', label="trajectoire de la Terre") #trajectoire de la Terre
    plt.plot(x, y, '-r', label="trajectoire de la Lune")
    plt.title("Trajectoire de la Terre et de la Lune")
    plt.legend()
    #plt.savefig('figure.jpg')
    plt.show()


### Programme principal

G=6.67259e-11 # SI
M_S=1.989e30 # kg
M_T=5.972e24 # kg
M_L=7.3477e22 # kg
dS_T=149.597870700e9 # m
dT_L=381.5e6 # m
v0=np.sqrt(G*M_S/dS_T) # = 29785.294991204377 m/s , vitesse pour la Terre
v0L=1022 # m/s , vitesse de la Lune

dt=100 # s
N=315582 # *dt -> nombre de secondes dans une année (une révolution de la Terre), 365.256363051 est le nombre de jours dans une année.

tE, X_T, Y_T, X_L, Y_L = Euler(fT_Lx, fT_Ly, fS_Tx ,fS_Ty, dS_T, dS_T+dT_L, v0, v0+v0L, dt, N)


tmp_2=time.time()

tmp=tmp_2-tmp_1
print('\n ', tmp, "secondes pour exécuter le programme")



plot_traj(X_L, Y_L, X_T, Y_T, 2, "une")

plt.plot(tE, np.sqrt(X_L**2+Y_L**2) - np.sqrt(X_T**2+Y_T**2) )
plt.title("Distance Terre-Lune en fonction du temps.")
plt.xlabel('temps (en s)')
plt.ylabel('distance (en m)')
plt.show()


