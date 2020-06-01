##### Olivier Roth 
##### SN4 
##### Equation trajectoire M



### Importation des modules

import numpy as np
import matplotlib.pyplot as plt
import progressbar # pour une barre de chargement du programme (à utiliser si le module a été téléchargé)
import time

tmp_1=time.time()



### Constantes

G=6.67259e-11 # SI
M_S=1.989e30 # kg
M_T=5.972e24 # kg
dS_T=149.597870700e9 # m
d1=1.4917449889577318e9 # m, distance Terre - point L1
d2=-1.501728936924601e9 # m, distance Terre - point L2
m= 100  # kg, masse de l'objet
v0=np.sqrt(G*M_S/dS_T) # vitesse de la Terre
v0m1=np.sqrt(G*(M_S+M_T)/(dS_T+d1)) # vitesse au point L1
v0m2=np.sqrt(G*(M_S+M_T)/(dS_T+d2)) # vitesse au point L2



### Fonctions

fx = lambda x, y , X_T, Y_T : -G*m*M_S*x/((x**2+y**2)**(3/2)) - G*m*M_T*(x-X_T)/(((x-X_T)**2+(y-Y_T)**2)**(3/2)) # force en x exercée sur un objet de masse m
fy = lambda x, y, X_T, Y_T : -G*m*M_S*y/((x**2+y**2)**(3/2)) - G*m*M_T*(y-Y_T)/(((x-X_T)**2+(y-Y_T)**2)**(3/2)) # force en y exercée sur un objet de masse m

fS_Tx= lambda x_t, y_t : -G*M_T*M_S*x_t/((x_t**2+y_t**2)**(3/2)) # Force du Soleil sur la Terre en x
fS_Ty= lambda x_t, y_t : -G*M_T*M_S*y_t/((x_t**2+y_t**2)**(3/2)) # Force du Soleil sur la Terre en y


def Euler(fx, fy, fS_Tx, fS_Ty, dS_T, d0, v0, v1, dt, N) : # fonction principale
    ttab=np.zeros(N)
    xtab=np.zeros(N)
    x2tab=np.zeros(N)
    ytab=np.zeros(N)
    y2tab=np.zeros(N)
    vxtab=np.zeros(N)
    vx2tab=np.zeros(N)
    vytab=np.zeros(N)
    vy2tab=np.zeros(N)
    
    xtab[0]=dS_T # conditions initiales
    x2tab[0]=dS_T-d0
    vytab[0]=v0
    vy2tab[0]=v1
    
    bar = progressbar.ProgressBar(widgets=[progressbar.Percentage(),progressbar.Bar(),], max_value=N).start() # barre de chargement
    
    for i in range(N-1) :
        ttab[i+1] = ttab[i]+dt
        xtab[i+1] = xtab[i]+dt*vxtab[i]
        x2tab[i+1] = x2tab[i]+dt*vx2tab[i]
        ytab[i+1] = ytab[i]+dt*vytab[i]
        y2tab[i+1] = y2tab[i]+dt*vy2tab[i]
        vxtab[i+1] = vxtab[i]+dt*fS_Tx(xtab[i],ytab[i])/M_T
        vx2tab[i+1] = vx2tab[i]+dt*(fx(x2tab[i],y2tab[i], xtab[i], ytab[i]))/m
        vytab[i+1] = vytab[i]+dt*fS_Ty(xtab[i],ytab[i])/M_T
        vy2tab[i+1] = vy2tab[i]+dt*(fy(x2tab[i],y2tab[i], xtab[i], ytab[i]))/m
        bar += 1 # enlever le # si le module a été importé
    bar.finish() # enlever le # si le module a été importé
    return ttab, xtab, ytab, x2tab, y2tab



def v_optimale(n): # pour trouver la vitesse optimale
    tE, X_T, Y_T, X_M, Y_M = Euler(fx, fy, fS_Tx, fS_Ty, dS_T, d1, v0, v0m1, dt, N)
    a=np.zeros(n)
    b=10
    for i in range(n-1):
        b=100/10**(i)
        print('\n étape', i+1, 'sur', n-1) 
        print(' a=', a[i], 'b=', b)
        while X_M[-1]-X_M[0] >0 :
            tE, X_T, Y_T, X_M, Y_M = Euler(fx, fy, fS_Tx, fS_Ty, dS_T, d1, v0, v0m1-a[i], dt, N)
            print(' a=', a[i], '\n ')
            a[i]+=b
        tE, X_T, Y_T, X_M, Y_M = Euler(fx, fy, fS_Tx, fS_Ty, dS_T, d1, v0, v0m1-a[i]+2*b, dt, N)
        a[i+1]=a[i]-2*b
        print('\nChangement de pas')
    print('v=', v0m1-a[n-1], 'a=', a[n-1])
    return v0m1-a[n-1]


def instab(n, L, vl, v, a, b): # pour prouver l'instabilité de L1, L2
    for i in range(n):
        tE, X_T, Y_T, X_M, Y_M = Euler(fx, fy, fS_Tx, fS_Ty, dS_T, L+a*10**(i+1), v0, vl+v, dt, N)
        d=a*(np.sqrt(X_T**2+Y_T**2)-np.sqrt(X_M**2+Y_M**2))
        plt.plot(tE, d, label=b+'+'+str(10**(i+1))+' m')
    plt.xlabel('temps (s)')
    plt.ylabel("Ecart entre la Terre et l'objet en fonction de la distance initiale (m).")
    #plt.title("Graphiques de l'évolution de la distance Terre-objet en fct du temps pour plusieurs distances initiales.")
    plt.legend()
    plt.show()



def plot_traj(x, y, xt, yt, a, b) : # plot des résultats des équations différentielles
    plt.axis([-a*1e11,a*1e11,-a*1e11,a*1e11])
    plt.axis('equal')
    plt.plot(xt, yt, '-b', label="trajectoire de la Terre") #trajectoire de la Terre
    plt.plot(x, y, '-r', label="trajectoire de l'objet en L"+b)
    plt.title("Trajectoire de la Terre et d'un objet au point L"+b)
    plt.legend()
    #plt.savefig('figure.jpg')
    plt.show()


def plot_2traj(x, y, x2, y2, xt, yt, a):
    plt.axis([-a*1e11,a*1e11,-a*1e11,a*1e11])
    plt.axis('equal')
    plt.plot(xt, yt, '-b', label="trajectoire de la Terre") #trajectoire de la Terre
    plt.plot(x, y, '-r', label="trajectoire de l'objet en L1")
    plt.plot(x2, y2, '-g', label="trajectoire de l'objet en L2")
    plt.title("Trajectoire de deux objets ")
    plt.legend()
    plt.show()


def plot_anim(N, a, x, y) : # pour un plot animé
    plt.ion()
    n=0
    while n<N:
	    #plt.clf() # pour 'clear' les points précédents
	    plt.axis([-a*1e11,a*1e11,-a*1e11,a*1e11])
	    plt.plot(x[n],y[n],'xr')
	    plt.draw()
	    plt.pause(1e-2)
	    n+=10000
    plt.ioff()
    plt.show()	






### Programme principal

dt=100 # s
N=315582 # *dt -> nombre de secondes dans une année (une révolution de la Terre), avec 365.256363051 le nombre de jours dans une année.


#vL1=v_optimale(10) # pour trouver la vitesse optimale en L1, à 10**(-6) près
#tE, X_T, Y_T, X_M, Y_M = Euler(fx, fy, fS_Tx, fS_Ty, dS_T, d1, v0, vL1, dt, N) # orbite en L1

v1=-149.175379 # correction de vitesse pour L1
v2=147.84918 # correction de vitesse pour L2


tE, X_T, Y_T, X_M, Y_M = Euler(fx, fy, fS_Tx, fS_Ty, dS_T, d1, v0, v0m1+v1, dt, N) # orbite en L1
#tE, X_T, Y_T, X_M2, Y_M2 = Euler(fx, fy, fS_Tx, fS_Ty, dS_T, d2, v0, v0m2+v2, dt, N) # orbite en L2



tmp_2=time.time()

tmp=tmp_2-tmp_1
print('\n ', tmp, "secondes pour exécuter le programme")


#instab(3, d1, v0m1, v1, 1, 'L1') # pour prouver l'instabilité de L1, (fait un plot)
#instab(3, d2, v0m2, v2, -1, 'L2') # pour prouver l'instabilité de L2, (fait un plot)



### Plot

#plot_traj(X_M, Y_M, X_T, Y_T, 2, "1") # objet en L1
#plot_traj(X_M2, Y_M2, X_T, Y_T, 2, "2") # objet en L2

#plot_2traj(X_M, Y_M, X_M2, Y_M2, X_T, Y_T, 2) # traj en L1 et L2


plot_anim(N, 2, X_M, Y_M) # pour un graph animé


    




