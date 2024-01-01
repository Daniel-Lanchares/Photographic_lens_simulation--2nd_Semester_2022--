# -*- coding: utf-8 -*-
"""
Created on Sat Feb 19 11:47:44 2022

@author: danie
"""

#Trabajo Óptica V3

import numpy as np
import matplotlib.pyplot as plt


class Intercara:
    def __init__(self, n0, n1, R01, config={}): #Ojo, Si lente Cóncava (Divergente) R01 = -|R01|
        
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
        
        
        self.config_predet()
        
        for attr, val in config.items():
            setattr(self, attr, val)
    
    def config_predet(self):
        self.pos = 0
        self.h = 5
    
    def dibuja(self, ejes):
        if self.R01 == np.inf: ejes.plot(self.pos*np.ones(20), np.linspace(-self.h, self.h, 20), 'b')
        else:
            if self.h/abs(self.R01) >1: 
                ang_max = np.pi/2
                plt.plot(self.pos*np.ones(2), [abs(self.R01),self.h], 'b-')
                plt.plot(self.pos*np.ones(2), [-abs(self.R01),-self.h], 'b-')
            else: ang_max=np.arcsin(self.h/abs(self.R01))
            ang = np.linspace(-ang_max, ang_max, 100)
            ejes.plot(self.R01*np.cos(ang) + (self.pos-self.R01*np.cos(ang_max)), self.R01*np.sin(ang), 'b-')
    
    def dibuja_planos(self, ejes):
        
        Ppos = ejes.plot((self.pos)*np.ones(20), np.linspace(-self.R01, self.R01, 20), 'b-')
        if self.fi == np.inf: pass
        else: 
            PFi = ejes.plot((self.fi+self.pos)*np.ones(20), np.linspace(-self.R01, self.R01, 20), 'g-')
            PF0 = ejes.plot((self.pos-self.f0)*np.ones(20), np.linspace(-self.R01, self.R01, 20), 'k-')
        

class Haz:
    def __init__(self, origen, pendiente, config={}): #TERMINAR y integrar dibuja rayo
        
        self.haz = []
        self.haz_dibujo = []#No incluye los haces de los traslados
        self.posiciones = [origen[0],]
        if isinstance(pendiente, (int, float)): 
            self.n =  len(origen[1])
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

def sistema(Haces):
    for Haz in Haces:
        Haz.trasl(I1.pos)
        for i in range( len(sist)):
            Haz.pasa_obj(sist[i])
        
            if i==len(sist)-1: 
                Haz.trasl(46.4)
            else: 
                Haz.trasl(sist[i+1].pos-sist[i].pos)


def dibuja_sist(ejes):
    for I in sist:
        I.dibuja(ejes)

def simulacion(sist, haces, ejes):
    Mf = np.identity(2)
    for i in range(len(sist)):
        Mf = sist[i].M @ Mf
        if i == len(sist)-1:
            continue
        else:
            d = sist[i+1].pos-sist[i].pos
            D = np.array([[1,d],[0,1]])
            Mf = D @ Mf
    
    fi = -1/Mf[1][0]
    di = Mf[0][0]*fi
    
    do = -Mf[1][1]/Mf[1][0]
    fo = Mf[0][0]*do + Mf[0][1]
    '''
    print(Mf)
    print(f'distancia focal imagen (fi): {fi}')
    print(f'distancia focal imagen (di): {di}')
    print(f'Punto de Convergencia: {di + sist[-1].pos}')
    print(f'distancia focal objeto (fo): {fo}')
    print(f'distancia focal objeto (do): {do}')
    '''
    
    
    

    ejes.set_aspect('equal', 'box')
    
    #Planos principales
    y_planos= np.linspace(-35/2, 35/2, 20)
    
    ejes.plot((di + sist[-1].pos)*np.ones(20), y_planos, 'r', label='Plano Focal Imagen')
    ejes.plot((di + sist[-1].pos - fi)*np.ones(20), y_planos, 'k', label='Plano Principal')
    ejes.plot((-do + sist[-1].pos)*np.ones(20), y_planos, 'g', label='Plano Focal Objeto')
    
    ejes.plot((pos_sensor)*np.ones(20),y_planos , 'k', label='Sensor')
    ejes.plot(np.linspace((di + sist[-1].pos),(pos_sensor), 20), -35**np.ones(20), 'k' )
    
    
    sistema(haces)
    dibuja_sist(ejes)
    ejes.text( 150, -40,f'distancia focal: {round(pos_sensor -(di + sist[-1].pos), 2)}')
    
    
    
    for H in haces: H.dibuja_haz(ejes)
    
    if rayos_prolongados:
        p_p = H2.haz[-1][-1]
        pos= H2.posiciones[-1]
        x = np.arange(10, 60)
        y = p_p[0] + p_p[1]*(x-pos)
        plt.plot(np.arange(60), H2.haz[0][-1][0]*np.ones(60))
        plt.plot(x, y)
    
    ejes.plot(np.linspace(0, max([haz.posiciones[-1] for haz in haces]), 50), np.zeros(50), 'k-')
    
    #ejes.set_xlim(-lx/2,lx/2)
    ejes.set_ylim(-18*factor,18*factor)
    ejes.set_aspect('equal', 'box')
    plt.legend()
    #plt.show()

'''
Construir  sistemas ópticos más complejos, 
y obtener sus matrices para conseguir puntos y distancias focales
'''

factor = 3
    

I1 = Intercara(n0=1,  n1=1.5,R01= -70, config={'h':50})

I2 = Intercara(n0=1,  n1=1.5,  R01= np.inf, config={'h':40})
I3 = Intercara(n0=1,  n1=1.5,R01= -7*factor, config={'h':7*factor}) 

I4 = Intercara(n0=1.5,n1=1,  R01= -12*factor, config={'h':7*factor}) 
I5 = Intercara(n0=1,  n1=1.5,R01= 12*factor, config={'h':7*factor})

I6 = Intercara(n0=1.5,n1=1,  R01= 6*factor, config={'h':5*factor})
I7 = Intercara(n0=1,  n1=1.5,R01= 10*factor, config={'h':5*factor}) 

I8 = Intercara(n0=1.5,n1=1,  R01= -15*factor, config={'h':6*factor}) #R01<0 Convergente, R01>0 Divergente
I9 = Intercara(n0=1.5,  n1=1,R01=np.inf, config={'h':6*factor})


sist=[I1, I2, I3, I4, I5, I6, I7, I8, I9]


#l=0
rayos_prolongados = False

fig = plt.figure(figsize=(16,8))
ejes = fig.add_subplot(1,1,1)

for l in [0,]:#np.linspace(0, 6, 100):
    
    I1.pos = 10*factor
    I2.pos = 15*factor#11
    
    I3.pos = 23*factor#16 + l  
    I4.pos = (35-1.2*l)*factor#20 + l 
    
    I5.pos = (38-1.2*l)*factor#25 + 1.2*l
    I6.pos = (45 - 1.3*l)*factor#30 + 1.2*l
    
    I7.pos = (49 - 1.3*l)*factor#42
    I8.pos = 54*factor#50
    
    I9.pos = 55*factor#55
    
    pos_sensor = I9.pos + 46.4 #Montura Nikon F
    
    
    H1 = Haz([0,2], np.linspace(-0.2,0,3), config={
        'color':'r'})
    
    H2 = Haz([0, np.linspace(-6, 6, 15)], 0, config={
        'color': 'y'})
    
    H3 = Haz([0,0], np.linspace(-0.1,0.1,2), config={
        'color':'y'})
    
    H4 = Haz([0,-4], np.linspace(0,0.2,3), config={
        'color':'g'})
    
    haces = [H2,]
    

    ejes.cla()
    simulacion(sist, haces, ejes)
    
    plt.pause(0.05)