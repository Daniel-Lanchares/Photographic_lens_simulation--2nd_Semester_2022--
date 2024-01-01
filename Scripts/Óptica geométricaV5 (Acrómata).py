# -*- coding: utf-8 -*-
"""
Created on Sun Feb 20 20:05:58 2022

@author: danie
"""

#Acrómata

import numpy as np
from copy import deepcopy
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

class Intercara:
    def __init__(self, n0, n1, R01, config={}): #Ojo, Si lente Cóncava (Divergente) R01 = -|R01|
        
        self.R01 = R01
        self.n0 = n0
        self.n1 = n1
        self.n01 = lambda l: n0(l)/n1(l)
        if self.R01 == np.inf: 
            self.M = lambda l: np.array([[1,0],[0,self.n01(l)]])
            self.fi = np.inf
        else: 
            self.M = lambda l: np.array([[1,0],[(1-self.n01(l))/R01,self.n01(l)]])
            self.fi = lambda l: -1/self.M(l)[1,0]
        self.f0 = lambda l: self.n01(l)*self.fi(l)
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
    
    def dibuja_planos(self, ejes): #ahora mismo no funcionaría
        
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
        self.l = 500e-6 #mm
    
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
            nuevo_rayo.append(Obj.M(self.l)@rayo)
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

def abbe(n):
    B = 486.1327e-6
    Y = 587.5618e-6
    R = 656.2816e-6
    return (n(Y)-1)/(n(B)-n(R))

def n_doblete(inter, f):
    '''
    

    Parametros
    ----------
    inter : TYPE Intercara
        Intercara a la que se añade la lente.
    f : TYPE float
        Focal de la lente a añadir

    Returns
    -------
    n : TYPE Float
        Índice de refracción ideal del doblete
    '''
    Y = 587.5618e-6
    V1 = abbe(inter.n0)
    V2 = -V1*inter.fi(Y)/f
    
    Vs = np.array([abs(float(Tabla[i][-1])-V2) for i in range(len(Tabla))])
    V_min = np.max(Vs)
    indice = np.where(Vs == (V_min))[0][0]
    n = Tabla[indice][1]
    return float(n)

def R_doblete(n1, n2, f, R1_1):
    Y = 587.5618e-6
    v1 = abbe(n1)
    v2 = abbe(n2)
    a = f*(v1-v2)
    R1_2 = 1/(v1/(a*(n1(Y)-1))-1/R1_1)
    R2_1 = -a*(n2(Y)-1)/v2
    return R1_2, R2_1


def sistema(Haces, sist):
    for Haz in Haces:
        Haz.trasl(sist[0].pos-Haz.posiciones[0])
        for i in range( len(sist)):
            Haz.pasa_obj(sist[i])
        
            if i==len(sist)-1: 
                Haz.trasl(d_final)
            else: 
                Haz.trasl(sist[i+1].pos-sist[i].pos)


def colormap(n1):
    if   1.0 < n1(np.inf) <= 1.3: return 'aqua'
    elif 1.3 < n1(np.inf) <= 1.4: return 'teal'
    elif 1.4 < n1(np.inf) <= 1.5: return 'deepskyblue'
    elif 1.5 < n1(np.inf) <= 1.6: return 'blue'
    elif 1.6 < n1(np.inf) <= 1.7: return 'navy'

def H_cmap(l):
    L=l*1e6
    if         L <= 400: return 'darkviolet'
    elif 400 < L <= 450: return 'blue'
    elif 450 < L <= 500: return 'cyan'
    elif 500 < L <= 550: return 'yellowgreen'
    elif 550 < L <= 600: return 'yellow'
    elif 600 < L <= 650: return 'darkorange'
    elif 650 < L <= 700: return 'red'
    elif 700 < L       : return 'darkred'

def dibuja_sist(ejes, sist):
    for i,I in enumerate(sist):
        I.dibuja(ejes)
        
        #Rellenos
        if i < len(sist)-1 and I.n1 !=1:
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
        Mf = sist[i].M(Y) @ Mf
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
    
    
    

    
    #Planos principales
    
    y_planos = np.linspace(-35/2, 35/2, 20)
    y_planos2= np.linspace(-19.5/2, 19.5/2, 20)
    
    ejes.plot((di + sist[-1].pos)*np.ones(20), y_planos, 'r', label='Plano Focal Imagen')
    ejes.plot((di + sist[-1].pos - fi)*np.ones(20), y_planos, 'k', label='Plano Principal')
    #ejes.plot((-do + sist[-1].pos)*np.ones(20), y_planos, 'g', label='Plano Focal Objeto')
    
    if ejes == axs[1,1] or ejes == axs[0,1]: ejes.set_xlim(di + sist[-1].pos-2, di + sist[-1].pos+1)
    else: ejes.legend()
    '''
    ejes.plot((pos_sensor)*np.ones(20),y_planos , 'k')
    ejes.plot((pos_sensor)*np.ones(20),y_planos2 , 'r', label='Sensor')
    ejes.plot(np.linspace((di + sist[-1].pos),(pos_sensor), 20), -35**np.ones(20), 'k' )
    '''
    
    sistema(haces, sist)
    dibuja_sist(ejes, sist)
    #ejes.text( 150, -40,f'distancia focal: {round(pos_sensor -(di + sist[-1].pos), 2)}')
    
    
    
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
    
    ejes.set_aspect('equal', 'box')
    
    #plt.show()

'''
Tabla Wikipedia

Material	                   A	B (μm2)
Fused silica	            1.4580	0.00354
Borosilicate glass BK7	    1.5046	0.00420
Hard crown glass K5	        1.5220	0.00459
Barium crown glass BaK4	    1.5690	0.00531
Barium flint glass BaF10	1.6700	0.00743
Dense flint glass SF10	    1.7280	0.01342


Líneas Fraunhofer
Longitudes Referencia (nm) 

Azul        486.1327
Amarillo    587.5618
Rojo        656.2816
'''
Tabla = np.array([
    ['Crown Glass'          , 1.523, 58.6],
    ['High Index Glass'     , 1.600, 42],
    ['High Index Glass'     , 1.700, 39],
    ['Plastic CR-39'        , 1.590, 58],
    ['Flint extra-denso'    , 1.720, 29.3]])
B = 486.1327e-6
Y = 587.5618e-6
R = 656.2816e-6

l=0
d_final = 250

n_sil =     lambda l: 1.4580 +  3.54e-9/l**2
n_BK_7 =    lambda l: 1.5046 +  4.20e-9/l**2
n_K5 =      lambda l: 1.5220 +  4.59e-9/l**2
n_SF10 =    lambda l: 1.7280 + 13.42e-9/l**2
n_aire =    lambda l: 1 + 0/l**2


n_int = n_K5

sist = [#Lente sola
        Intercara(n0=n_aire,n1=n_int,  R01= -50, config={'h':15, 'pos':37}),
        Intercara(n0=n_int,  n1=n_aire,R01= 50, config={'h':15, 'pos':38}),
        #Intercara(n0=n_vid,n1=n_aire,  R01= np.inf, config={'h':15, 'pos':93}),
        
        ]

#n_corregido = n_doblete(sist[-1], -1.5*sist[-1].fi(Y))

print('K5')
print(f'Número de Abbe: {abbe(n_K5)}', f'n en d: {n_sil(Y)}')
print('SF10')
print(f'Número de Abbe: {abbe(n_SF10)}', f'n en d: {n_SF10(Y)}')
#print(R_doblete(n_sil, n_BK_7, 60, 50))

n_corr = n_SF10#lambda l: n_corregido

def cromatico(R, n2, Inter1, Inter2, l):
    phi1 = -(Inter1.n1(l)-1)*(1/Inter1.R01-1/Inter2.R01)
    V1 = abbe(Inter1.n1)
    phi2 = -(n2(l)-1)*(1/Inter2.R01-1/R)
    V2 = abbe(n2)
    return V1/phi1 + V2/phi2

sol = fsolve(cromatico, x0=105.83, args=(n_corr, sist[0], sist[1], Y))[0]

sist2 = [#Doblete
        Intercara(n0=n_aire,n1=n_int,  R01= -50, config={'h':15, 'pos':37}),
        Intercara(n0=n_int,  n1=n_corr, R01= 50, config={'h':15, 'pos':38}),
        Intercara(n0=n_corr,n1=n_aire,  R01= sol, config={'h':15, 'pos':43}), #R = 105.836
        
        ]



print(f'Socución: {sol}')
print(f'croma: {cromatico(sol, n_corr, sist[0],sist[1],Y)}')
print(1/(1/sist[0].fi(Y) + 1/sist[1].fi(Y)))


rayos_prolongados = False

fig, axs = plt.subplots(2,2,figsize=(16,9), gridspec_kw={'width_ratios': [7, 4]})

print(axs)

ejes = axs[0,0]
ejes2 = axs[1,0]

zoom1 = axs[0,1]
zoom2 = axs[1,1]

ejes.set_ylim(-20,20)
ejes.set_xlim(20, 150)

ejes2.set_ylim(-20,20)
ejes2.set_xlim(20, 150)

zoom1.set_ylim(-0.75,0.75)
zoom2.set_ylim(-0.75,0.75)

ejes.set_title(
    'Lente convergente de K5 (Crown)      $V_{K5}$ = %4.2f \n     $n_{K5}(\lambda) = 1.5220 + 4.59·10^{-9}/\lambda^2$ ' %abbe(n_K5))
ejes2.set_title(
    'Doblete con lente de SF10 (Flint)         $V_{SF10}$ = %4.2f \n    $n_{SF10}(\lambda) = 1.7280 + 13.42·10^{-9}/\lambda^2$ '%abbe(n_SF10))

zoom1.set_title('Zoom sobre plano focal imagen')
cr=format(cromatico(sol, n_corr, sist[0],sist[1],Y), '.2e')
zoom2.set_title(f'Cromático del doblete: {cr} mm')

for l in [0,]:#np.linspace(0, 6, 100):
    
    
    
    H1 = Haz([0, np.linspace(-9, 9, 3)], 0, config={
         'l': 350e-6})
    
    H2 =  Haz([0, np.linspace(-9, 9, 3)], 0, config={
         'l': 440e-6})
    
    H3 =  Haz([0, np.linspace(-9, 9, 3)], 0, config={
         'l': 540e-6})
    
    H4 =  Haz([0, np.linspace(-9, 9, 3)], 0, config={
         'l': 610e-6})
    
    H5 = Haz([0, np.linspace(-9, 9, 3)], 0, config={
        'l': 660e-6})
    Href = Haz([0, np.linspace(-9, 9, 3)], 0, config={
        'l': Y})
    
    haces = [H1, H2, H3, Href, H4, H5]
    for haz in haces: haz.color=H_cmap(haz.l)
    haces2 = [deepcopy(H) for H in haces]
    

    
    simulacion(sist, haces, ejes)
    simulacion(sist, haces, zoom1)
    simulacion(sist2, haces2, ejes2)
    simulacion(sist2, haces2, zoom2)
    
    
    plt.pause(0.05)