# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 18:48:27 2022

@author: danie
"""

#Trabajo Óptica

import numpy as np
import matplotlib.pyplot as plt


class Intercara:
    def __init__(self, n0, n1, R01, pos = 0): #Ojo, Si lente Cóncava (Divergente) R01 = -|R01|
        self.pos = pos
        self.R01 = R01
        self.n01 = n0/n1
        if self.R01 == np.inf: 
            self.M = np.array([[1,0],[0,self.n01]])
            self.fi = np.inf
        else: 
            self.M = np.array([[1,0],[(1-self.n01)/R01,self.n01]])
            self.fi = -1/self.M[1,0]
        self.f0 = self.n01*self.fi
        self.obj=[]
    
    def dibuja(self, ejes):
        if self.R01 == np.inf: ejes.plot(self.pos**np.ones(20), np.linspace(-5, 5, 20))
        else:
            ang = np.linspace(-np.pi/2, np.pi/2, 100)
            ejes.plot(self.R01*np.cos(ang) + self.pos, self.R01*np.sin(ang), 'b-')
    
    def dibuja_planos(self, ejes):
        
        Ppos = ejes.plot((self.pos)*np.ones(20), np.linspace(-self.R01, self.R01, 20), 'b-')
        if self.fi == np.inf: pass
        else: 
            PFi = ejes.plot((self.fi+self.pos)*np.ones(20), np.linspace(-self.R01, self.R01, 20), 'g-')
            PF0 = ejes.plot((self.pos-self.f0)*np.ones(20), np.linspace(-self.R01, self.R01, 20), 'k-')
        

class Haz:
    def __init__(self, origen, pendiente, n=11, config={}): #TERMINAR y integrar dibuja rayo
        
        self.haz = []
        self.haz_dibujo = []#No incluye los haces de los traslados
        self.posiciones = [origen[0],]
        if isinstance(pendiente, (int, float)): 
            self.n =  n
            pendiente = pendiente*np.ones(self.n)
        else: 
           self.n = len(pendiente)
           
        
        if isinstance(origen[1], (int, float)): y = origen[1]*np.ones(self.n)
        else: y = origen[1]
        
        self.rayos = np.array(np.column_stack((y, pendiente)))
        self.haz.append(self.rayos) #Cada vez que se multiplique por una matriz, que añada los transformados para luego pintar.
        self.haz_dibujo.append(self.rayos)
        self.config_predet()
        
        for attr, val in config.items():
            setattr(self, attr, val)
    def config_predet(self):
        self.color='r'
    
    def trasl(self, d):
        M = np.array([[1,d],[0,1]])
        nuevo_rayo = []
        for rayo in self.haz[-1]:
            nuevo_rayo.append(M@rayo)
        self.haz.append(np.array(nuevo_rayo))
        self.posiciones.append(self.posiciones[-1]+d)

    def pasa_obj(self, Obj):
        nuevo_rayo = []
        for rayo in self.haz[-1]:
            nuevo_rayo.append(Obj.M@rayo)
        self.haz.append(np.array(nuevo_rayo))
        self.haz_dibujo.append(np.array(nuevo_rayo))

    def dibuja_rayos(self, rayos, interv, ejes):
        color = self.color
        for rayo in rayos:
             #Dibuja un rayo recto en un intervalo de x dado
             x = np.linspace(interv[0],interv[1],100)
             ejes.plot(x, rayo[1]*(x-interv[0])+rayo[0], color)
    
    def dibuja_haz(self, ejes):
        for i in range(len(self.haz_dibujo)):
            self.dibuja_rayos(self.haz_dibujo[i], self.posiciones[i:i+2], ejes)

'''
class Lente_d:
    def __init__(self, R, d, n, Convergente=True, config={}):
        self.R = R
        self.d = d
        self.conv = Convergente
        self.n = n
        
        self.config_predet()
        
        for attr, val in config.items():
            setattr(self, attr, val)
        
        self.n0
        
        if self.conv: R = -R
        
        self.inter_a = Intercara(n0=self.n0, n1=self.n1, R01 = R)
        self.inter_b = Intercara(n0=self.n0, n1=self.n1, R01 = -R)
        self.obj=[self.inter_a, self.inter_b]
    def config_predet(self):
        self.n0 = 1
    def actualiza(self, pos): #Mueve todos los subobjetos a la vez
        self.pos= pos
        self.inter_a.pos = pos-self.d/2
        self.inter_b.pos = pos+self.d/2
    
    def dibuja(self,ejes):
        self.inter_a.dibuja(ejes)
        self.inter_a.dibuja_planos(ejes)
        self.inter_b.dibuja(ejes)
        self.inter_b.dibuja_planos(ejes)
'''
'''
Construir  sistemas ópticos más complejos, 
y obtener sus matrices para conseguir puntos y distancias focales
'''


Lente1 = Intercara(n0=1,  n1=1.5,R01=-5)
Lente2 = Intercara(n0=1.5,n1=1,  R01= 5) #R01<0 Convergente, R01>0 Divergente
Lente3 = Intercara(n0=1,  n1=1.5,R01= 5)
l = 15 #Lente Delgada si l = Lente1.pos??


Lente1.pos = 10
Lente2.pos = l
Lente3.pos = l+15
def sistema(Haces):
    for Haz in Haces:
        Haz.trasl(Lente1.pos)
    
        Haz.pasa_obj(Lente1)
    
        Haz.trasl(Lente2.pos-Lente1.pos)
    
        Haz.pasa_obj(Lente2)
    
        Haz.trasl(Lente3.pos-Lente2.pos)
        
        Haz.pasa_obj(Lente3)
        Haz.trasl(5)


H1 = Haz([0,2], np.linspace(-0.3,0.3,11), config={
    'color':'r'})

H2 = Haz([0, np.linspace(-2, 2, 11)], 0, config={
    'color': 'y'})


fig = plt.figure(figsize=(16,8))
ejes = fig.add_subplot(1,1,1)
ejes.set_aspect('equal', 'box')


Lente1.dibuja(ejes)
Lente1.dibuja_planos(ejes)

Lente2.dibuja(ejes)
Lente2.dibuja_planos(ejes)

Lente3.dibuja(ejes)
sistema([H1,H2])

H1.dibuja_haz(ejes)
#H2.dibuja_haz(ejes)

ejes.plot(np.linspace(0, 45, 50), np.zeros(50), 'k-')
plt.show()
