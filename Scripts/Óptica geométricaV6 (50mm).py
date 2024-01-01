# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 12:10:17 2022

@author: danie
"""


import numpy as np
import matplotlib.pyplot as plt


class Intercara:
    def __init__(self, n0, n1, R01, config={}): #Ojo, Si lente Cóncava (Divergente) R01 = -|R01|
        
        self.R01 = R01
        self.n0 = n0
        self.n1 = n1
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
        
        if self.R01 != np.inf: self.corte = self.pos+self.R01*(1-np.cos(np.arcsin(self.h/self.R01))) #Donde corta con el eje óptico
        else: self.corte = self.pos
    
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
            self.curva = [self.R01*np.cos(ang) + (self.pos-self.R01*np.cos(ang_max)) , self.R01*np.sin(ang)]
            ejes.plot(self.curva[0],self.curva[1] , 'b-')
    
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
        self.invertido = False
    
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

def sistema(Haces, sist):
    for Haz in Haces:
        if not Haz.invertido:
            Haz.trasl(sist[0].pos-Haz.posiciones[0])
            for i in range( len(sist)):
                Haz.pasa_obj(sist[i])
            
                if i==len(sist)-1: 
                    Haz.trasl(80)
                else: 
                    Haz.trasl(sist[i+1].pos-sist[i].pos)
        else: #No funciona
            Haz.trasl(sist[-1].pos-Haz.posiciones[0])
            for i in np.flip(np.arange(len(sist))):
                Haz.pasa_obj(sist[i])
            
                if i==len(sist)-1: 
                    Haz.trasl(80)
                else: 
                    Haz.trasl(sist[i+1].pos-sist[i].pos)


def colormap(n1):
    if   1.0 < n1 <= 1.3: return 'aqua'
    elif 1.3 < n1 <= 1.4: return 'teal'
    elif 1.4 < n1 <= 1.5: return 'deepskyblue'
    elif 1.5 < n1 <= 1.6: return 'blue'
    elif 1.6 < n1 <= 1.7: return 'navy'

def dibuja_sist(ejes, sist):
    for i,I in enumerate(sist):
        I.dibuja(ejes)
        
        #Rellenos
        if i <= len(sist)-1 and I.n1 !=1:
            I2 = sist[i+1]
            color = colormap(I.n1)
            if I.R01 < 0 and I2.R01 <0:
                x1=np.linspace(I.corte, min(I2.corte, I.pos), 50)
                x2=np.linspace(min(I2.corte, I.pos), max(I2.corte, I.pos), 50)
                x3=np.linspace(max(I2.corte, I.pos), I2.pos, 50)
                
                y1 = lambda x: np.sqrt(I.R01**2 -(x-(-I.R01 +I.corte))**2)
                y2 = lambda x: np.sqrt(I2.R01**2-(x-(-I2.R01+I2.corte))**2)
                y3 = max(I.h,I2.h)*np.ones_like(y2)
                
                ejes.fill_between(x1, y1(x1), color=color)
                if I2.corte < I.pos: ejes.fill_between(x2, y1(x2), y2(x2), color=color)
                else: ejes.fill_between(x2, y3, color=color)
                ejes.fill_between(x3, y2(x3), y3, color=color)
                
                ejes.fill_between(x1, -y1(x1), color=color)
                if I2.corte < I.pos: ejes.fill_between(x2, -y1(x2), -y2(x2), color=color)
                else: ejes.fill_between(x2, -y3, color=color)
                ejes.fill_between(x3, -y2(x3), -y3, color=color)
           
            elif I.R01*I2.R01 < 0:
                x1=np.linspace(min(I.corte, I.pos), max(I.corte, I.pos), 50)
                x2=np.linspace(max(I.corte, I.pos), min(I2.corte, I2.pos), 50)
                x3=np.linspace(min(I2.corte, I2.pos), max(I2.corte, I2.pos), 50)
                
                y1 = lambda x: np.sqrt(I.R01**2 -(x-(-I.R01 +I.corte))**2)#Sale un número inválido pero no da problemas
                y2 = max(I.h,I2.h)*np.ones_like(y1)
                y3 = lambda x: np.sqrt(I2.R01**2-(x-(-I2.R01+I2.corte))**2)
                
                if I.R01 <0:
                    ejes.fill_between(x1, y1(x1), color=color)
                    ejes.fill_between(x2, y2, color=color)
                    ejes.fill_between(x3, y3(x3), color=color)
                    
                    ejes.fill_between(x1, -y1(x1), color=color)
                    ejes.fill_between(x2, -y2, color=color)
                    ejes.fill_between(x3, -y3(x3), color=color)
                else:
                    ejes.fill_between(x1, y1(x1), y2, color=color)
                    ejes.fill_between(x2, y2, color=color)
                    ejes.fill_between(x3, y3(x3), y2, color=color)
                    
                    ejes.fill_between(x1, -y1(x1), -y2, color=color)
                    ejes.fill_between(x2, -y2, color=color)
                    ejes.fill_between(x3, -y3(x3), -y2, color=color)
            elif I.R01 > 0 and I2.R01 > 0:
                x1=np.linspace(I.pos, min(I2.pos, I.corte), 50)
                x2=np.linspace(min(I2.pos, I.corte), max(I2.pos, I.corte), 50)
                x3=np.linspace(max(I2.pos, I.corte), I2.corte, 50)
                
                y1 = lambda x: np.sqrt(I.R01**2 -(x-(-I.R01 +I.corte))**2)
                y2 = lambda x: np.sqrt(I2.R01**2-(x-(-I2.R01+I2.corte))**2)
                y3 = max(I.h,I2.h)*np.ones_like(y2)
                
                ejes.fill_between(x1, y1(x1), y3, color=color)
                if I2.pos < I.corte: ejes.fill_between(x2, y1(x2), y2(x2), color=color)
                else: ejes.fill_between(x2, y3, color=color)
                ejes.fill_between(x3, y2(x3), color=color)
                
                ejes.fill_between(x1, -y1(x1), -y3, color=color)
                if I2.pos < I.corte: ejes.fill_between(x2, -y1(x2), -y2(x2), color=color)
                else: ejes.fill_between(x2, -y3, color=color)
                ejes.fill_between(x3, -y2(x3), color=color)
    
    
    

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
    
    print(Mf)
    print(f'distancia focal imagen (fi): {fi}')
    print(f'distancia focal imagen (di): {di}')
    print(f'Punto de Convergencia: {di + sist[-1].pos}')
    print(f'distancia focal objeto (fo): {fo}')
    print(f'distancia focal objeto (do): {do}')
    
    
    
    

    ejes.set_aspect('equal', 'box')
    
    #Planos principales
    y_planos = np.linspace(-35/2, 35/2, 20)
    y_planos2= np.linspace(-19.5/2, 19.5/2, 20)
    
    ejes.plot((di + sist[-1].pos)*np.ones(20), y_planos, 'r', label='Plano Focal Imagen')
    ejes.plot((di + sist[-1].pos - fi)*np.ones(20), y_planos, 'k', label='Plano Principal')
    ejes.plot((-do + sist[1].pos)*np.ones(20), y_planos, 'g', label='Plano Focal Objeto')
    
    ejes.plot((pos_sensor)*np.ones(20),y_planos , 'k')
    ejes.plot((pos_sensor)*np.ones(20),y_planos2 , 'r', label='Sensor')
    ejes.plot(np.linspace((di + sist[-1].pos),(pos_sensor), 20), -35**np.ones(20), 'k' )
    
    
    sistema(haces, sist)
    dibuja_sist(ejes, sist)
    ejes.text( 150, -40,f'distancia focal: {round(pos_sensor -(di + sist[-1].pos), 2)}')
    ejes.set_xlabel('Altura (mm)')
    ejes.set_ylabel('Longitud (mm)')
    
    
    for H in haces: H.dibuja_haz(ejes)
    
    if rayos_prolongados:
        p_p = H2.haz[-1][-1]
        pos= H2.posiciones[-1]
        x = np.arange(10, 60)
        y = p_p[0] + p_p[1]*(x-pos)
        plt.plot(np.arange(60), H2.haz[0][-1][0]*np.ones(60))
        plt.plot(x, y)
    
    #Eje Óptico
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
l=0
factor = 3
pos = 40

n1 = 1.713
n2 = 1.64831

sist = [#Grupo 1
        Intercara(n0=1,  n1=n1,R01= -88.702, config={'h':28, 'pos':pos + 6.3}),
        Intercara(n0=n1,  n1=1,R01= np.inf, config={'h':28, 'pos':pos+6.78}),
        
        Intercara(n0=1,  n1=n1,R01= -37.399, config={'h':22, 'pos':pos +6.78 + 7.50}),
        Intercara(n0=n1,  n1=1,R01= -61.587, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 }),
        
        
        Intercara(n0=1,  n1=n2,  R01= -78.525, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 + 0.4}),
        Intercara(n0=n2,  n1=1,R01= -31.585, config={'h':16, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 2}),
        
        #Grupo 2
        Intercara(n0=1,  n1=n2,R01= 35.134, config={'h':16, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 8}),
        Intercara(n0=n2,  n1=n1,R01= np.inf, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 8 + 6}),
        Intercara(n0=n1,  n1=1,R01= 46.274, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 8 + 6 + 1}),
        
        Intercara(n0=1,n1=n1,  R01= np.inf, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 8 + 6 + 1 + 5.5}),
        Intercara(n0=n1,  n1=1,R01= 80, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 8 + 6 + 1 + 5.5 + 2.62}),
        #Seguir por aquí
        
        ]




rayos_prolongados = False

fig = plt.figure(figsize=(16,8))
ejes = fig.add_subplot(1,1,1)

for l in [0,]:#np.linspace(0, 6, 100):
    
    
    
    pos_sensor = sist[-1].pos + 6 + 18 #Montura Sony E
    
    
    H1 = Haz([0,18], np.linspace(-0.6,0,3), config={
        'color':'tab:orange'})
    
    H2 = Haz([0, np.linspace(-20, 20, 9)], 0, config={
        'color': 'y'})
    
    H3 = Haz([0,0], np.linspace(-0.1,0.1,2), config={
        'color':'y'})
    
    H4 = Haz([-24.4658,-0], np.linspace(0,0.2,3), config={
        'color':'tab:purple'})
    
    haces = [H2,H1, H4] #Estudiar la distancia de enfoque
    

    ejes.cla()
    simulacion(sist, haces, ejes)
    
    plt.pause(0.05)
    print(H4.haz[-1])