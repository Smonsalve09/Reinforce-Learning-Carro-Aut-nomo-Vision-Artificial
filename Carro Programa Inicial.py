# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 14:47:13 2024

@author: 000010478
"""

import numpy as np
import matplotlib.pyplot as plt

class terreno:
    def __init__(self):
        da=10
        self.r=[]
        x=[10.25,9.56,9.68,10.89,12.46,14.03,17.78,19.35,18.75,16.45,12.82,10.4,5.69,1.81,-1.09,-6.05,-12.22,-14.64,-16.45,-18.87,-19.6,-19.23,-18.15,-14.4,-10.52,-7.26,-4.84,-2.06,-0.48,0.12,-0.73,-2.66,-5.69,-9.31,-11.61,-16.94,-18.87,-19.84,-19.72,-18.75,-16.57,-12.82,-9.07,-6.29,-2.78,0.12,2.42,5.08,8.47,10.65,12.82,16.21,19.48,20.69,22.02,22.86,22.98,21.05,19.72,17.66,15.6,14.27,12.7,11.01]
        y=[-1.0,5.1,9.81,13.54,15.49,15.97,17.6,20.68,25.23,28.31,27.82,27.66,27.18,27.01,27.01,27.18,27.5,27.18,26.85,25.06,21.49,17.76,13.54,11.75,11.59,12.73,14.84,13.54,10.13,6.56,2.66,-0.1,-0.58,-0.1,1.36,1.53,-0.1,-2.86,-6.92,-10.65,-11.95,-13.41,-14.55,-14.87,-14.06,-13.9,-13.9,-13.9,-16.98,-20.39,-21.85,-22.18,-21.04,-18.93,-15.03,-10.65,-7.4,-4.81,-3.18,-2.69,-3.34,-5.13,-4.32,-3.99]
        #for ang in np.arange(0,360,da):
        #    self.r.append(np.array([10*np.cos(ang*np.pi/180), 20*np.sin(ang*np.pi/180)]))
        for xi, yi in zip(x,y):
            self.r.append(np.array([xi,yi]))
        self.a=2 #Ancho de la via a cada lado del eje
        
    def tipo(self, ri): #Verifica si es terreno o camino
        dist=[np.linalg.norm(ri-r) for r in self.r]
        ind1=np.argmin(dist)
        if(dist[ind1]>self.a*2):
            return 1
        ind2=ind1+1 if ind1+1<len(self.r) else 0
        ind0=ind1-1 if ind1>0 else len(self.r)-1
        return max([self.incluye(ri,ind0), self.incluye(ri,ind1), self.incluye(ri,ind2)])
    def incluye(self, r, ind):
        #print('Metodo incluye: r=',r,', ind=',ind)
        th=6
        C1=self.r[ind]
        if(ind+1<len(self.r)):
            C2=self.r[ind+1]
        else:
            C2=self.r[0]
        #plt.figure()
        #plt.plot([C1[0],C2[0]],[C1[1],C2[1]],'o-k')
        #plt.plot(r[0], r[1], 'xr')
        #xc, yc=np.cos(np.arange(0,2*np.pi, 0.1)), np.sin(np.arange(0,2*np.pi, 0.1))
        #plt.plot(C1[0]+xc*self.a/th, C1[1]+yc*self.a/th, '-r')
        #plt.plot(C1[0]+xc*self.a, C1[1]+yc*self.a, '-b')
        #plt.plot(C2[0]+xc*self.a/th, C2[1]+yc*self.a/th, '-r')
        #plt.plot(C2[0]+xc*self.a, C2[1]+yc*self.a, '-b')
        #plt.axis('equal')
        gamma=C2-C1
        dist=np.linalg.norm(gamma)
        gamma=gamma/dist
        a=np.dot(r-C1,gamma)
        #print('a=',a, ', dist=',dist)
        if(a<0 or a>dist):
            if(np.linalg.norm(r-C1)<self.a/th or np.linalg.norm(r-C2)<self.a/th):
                #print('Cercanía interna a un extremo')
                return 3
            if(np.linalg.norm(r-C1)<self.a or np.linalg.norm(r-C2)<self.a):
                #print('Cercanía externa a un extremo')
                return 2
            #print('Está por fuera de la región de cobertura')
            return 1
        dp=np.linalg.norm(r-C1-a*gamma)
        #print('Distancia perpendicular: ',dp)
        #if(dp<=self.a/th):
        #    print(dp)
        if(dp<=self.a/th):
            return 3
        if(dp<=self.a):
            return 2
        return 1
    def dibujar(self):
        x=[ri[0] for ri in self.r]
        y=[ri[1] for ri in self.r]
        plt.plot(x,y,'o-k')
    def visor(self, robs=np.array([0,0]), mu=np.array([1,0]), D=2, H=2, W=3, dL=0.2, h=1):
        x,y=robs
        xp, yp=np.meshgrid(np.arange(-W,W+dL, dL), np.arange(-H,H+dL, dL))
        z=np.zeros_like(xp)
        Xs=np.zeros_like(xp)
        Ys=np.zeros_like(xp)
        for i in range(len(xp)):
            for j in range(len(xp[0])):
                if(yp[i,j]>=0):
                    z[i,j]=0
                else:
                    alfa=-h/yp[i,j]
                    xs=x+alfa*D*mu[0]+alfa*xp[i,j]*mu[1]
                    Xs[i,j]=xs
                    ys=y+alfa*D*mu[1]-alfa*xp[i,j]*mu[0]
                    Ys[i,j]=ys
                    z[i,j]=self.tipo(np.array([xs,ys]))
        return xp, yp, z, Xs, Ys

#Simulación del movimiento del carro. 
class carro:
    def __init__(self, dt=0.2, r=np.array([0,0]), fi=np.pi/2, theta=0, L=2, W=2, pista=None):
        self.dt=dt
        self.r=r #Posición del carro
        self.fi=fi #Dirección del movimiento del carro
        self.mu=np.array([np.cos(self.fi), np.sin(self.fi)])
        self.theta=theta #Ángulo de la dirección del carro. 
        self.L=L #Longitud del carro hasta la dirección
        self.W=W #Velocidad del carro en m/s
        self.pista=pista
    def cambio_theta(self, dtheta):
        self.theta=np.min([np.pi/2, np.max([-np.pi/2, self.theta+dtheta])])
    def mover(self):
        if(self.theta==0):
            self.r=self.r+self.W*self.mu*self.dt
        else:
            v=np.array([np.sin(self.fi), -np.cos(self.fi)])
            self.r=self.r+self.mu*self.L/np.sin(self.theta)*np.sin(self.W*self.dt*np.sin(self.theta)/self.L)+v*self.L/np.sin(self.theta)*(1-np.cos(self.W*self.dt*np.sin(self.theta)/self.L))
            self.fi=self.fi+self.dt*self.W*np.sin(self.theta)/self.L
            self.mu=np.array([np.cos(self.fi), np.sin(self.fi)])
    
    
    
colores={0:np.array([133, 193, 233]), 1:np.array([25, 111, 61]), 2: np.array([93, 109, 126]), 3: np.array([241, 196, 15])}
pista=terreno()
x,y=np.meshgrid(np.linspace(-30,30,121), np.linspace(-30,30,121))
z=np.zeros_like(x)
for i in range(len(x)):
    for j in range(len(x[0])):
        z[i,j]=pista.tipo(np.array([x[i,j], y[i,j]]))
plt.close('all')
zc=np.zeros((z.shape[0], z.shape[1], 3), dtype=int)
for i in range(len(zc)):
    for j in range(len(zc[0])):
        zc[i,j,:]=colores[z[i,j]]
plt.imshow(zc, origin='lower')    
plt.axis('equal')
plt.figure(figsize=(12,6))
cambiar=True
ind=-1
dS=0.5
robs=pista.r[0]
"""while(ind<len(pista.r)-1):
    if(cambiar):
        cambiar=False
        ind+=1
        C1=robs
        C2=pista.r[ind+1]
        dist=np.linalg.norm(C2-C1)
        mu=(C2-C1)/dist
    else:
        robs=robs+dS*mu
        dist-=dS
        if(dist<=0):
            cambiar=True
    D=1.5
    xp, yp, zp, xs, ys=pista.visor(robs=robs, mu=mu, D=D, dL=0.1)
    zc=np.zeros((zp.shape[0], zp.shape[1], 3), dtype=int)
    for i in range(len(zc)):
        for j in range(len(zc[0])):
            zc[i,j,:]=colores[zp[i,j]]
    plt.subplot(1,2,1)
    plt.gca().cla()
    plt.contourf(x,y,z)    
    plt.plot(robs[0], robs[1], 'xr')
    plt.plot([robs[0], robs[0]+mu[0]*D], [robs[1], robs[1]+mu[1]*D], '-r')
    plt.axis('equal')
    plt.subplot(1,2,2)
    plt.gca().cla()
    plt.imshow(zc, origin='lower')
    plt.pause(0.01)
"""
#Caso 2, carro moviendose según acciones
C=carro()
for k in range(200):
    D=1.5
    xp, yp, zp, xs, ys=pista.visor(robs=C.r, mu=C.mu, D=D, dL=0.025)
    zc=np.zeros((zp.shape[0], zp.shape[1], 3), dtype=int)
    for i in range(len(zc)):
        for j in range(len(zc[0])):
            zc[i,j,:]=colores[zp[i,j]]
    plt.subplot(1,2,1)
    plt.gca().cla()
    plt.contourf(x,y,z)    
    plt.plot(C.r[0], C.r[1], 'xr')
    plt.plot([C.r[0], C.r[0]+C.mu[0]*D], [C.r[1], C.r[1]+C.mu[1]*D], '-r')
    plt.axis('equal')
    plt.subplot(1,2,2)
    plt.gca().cla()
    plt.imshow(zc, origin='lower')
    plt.title(str(k))
    plt.pause(0.01)
    C.mover()
    """if(k>50 and k<70):
        C.cambio_theta(2*np.pi/180)
    if(k>80 and k<120):
        C.cambio_theta(-2*np.pi/180)
    if(k>120 and k<140):
        C.cambio_theta(2*np.pi/180)
    """
    if(k==50):
        C.cambio_theta(30*np.pi/180)
    if(k==80):
        C.cambio_theta(-30*np.pi/180)
    if(k==100):
        C.cambio_theta(60*np.pi/180)
             
        