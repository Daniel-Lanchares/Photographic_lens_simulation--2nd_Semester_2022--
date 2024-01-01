# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 09:18:58 2022

@author: danie
"""

#90 mm (enfoque))
import numpy as np
import matplotlib.pyplot as plt
from Simulador.Intercara import *
from Simulador.Diafragma import *
from Simulador.Haz import *





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
    
    
    

def simulacion(sist, haces, ejes, j):
    global pos
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
    ejes2.set_aspect('equal', 'box')
    
    #Planos principales
    y_planos = np.linspace(-25, 25, 20)
    y_planos2= np.linspace(-12, 12, 20)
    
    ejes.plot((di + sist[-1].pos)*np.ones(20), y_planos, 'r', label='Plano Focal Imagen')
    ejes.plot((di + sist[-1].pos - fi)*np.ones(20), y_planos, 'k--', label='Plano Principal')
    ejes.plot((-do + sist[1].pos)*np.ones(20), y_planos, 'g', label='Plano Focal Objeto')
    
     
    Planoi = di + sist[-1].pos
    for Obj in sist: 
        if type(Obj) == Diafragma: 
            Obj.pos = Planoi -fi
            fn = round(fi/(2*(Obj.h)), 1)
            ejes.text( round(di + sist[-1].pos)-30, -45,f'Número f: {fn}')
    Diff.append(Planoi)
    
    ejes.plot((pos_sensor)*np.ones(20),y_planos*0.8 , 'k')
    ejes.plot((pos_sensor)*np.ones(20),y_planos2 , 'tab:orange', label='Sensor')
    #ejes.plot(np.linspace((di + sist[-1].pos),(pos_sensor), 20), -35**np.ones(20), 'k' )
    ejes.plot(np.linspace((di + sist[-1].pos - fi),(di + sist[-1].pos), 20), -35**np.ones(20), 'k' )
    

    #Ángulo de visión general
    PP = di + sist[-1].pos - fi
    m = 12/(pos_sensor-PP)
    n = -m*PP
    x = np.linspace(sist[0].pos-30, pos_sensor + 30, 200)
    ejes.plot(x, m*x+n, 'r--')
    ejes.plot(x, -m*x-n, 'r--')
    
    sistema(haces, sist)
    dibuja_sist(ejes, sist)
    ejes.text( round(di + sist[-1].pos)-30, -40,f'distancia focal: {round(fi)}')
    #ejes.text( round(di + sist[-1].pos)-20, -50,f'Plano focal: {round(di + sist[-1].pos-pos, 2)}')#pos_sensor -(di + sist[-1].pos), 2)
    ejes.set_ylabel('Altura (mm)')
    ejes.set_xlabel('Longitud (mm)')
    ejes.set_title('Objetivo Fotográfico 90mm (Doble Gauss)')
    
    
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
    ejes.set_xlim(sist[0].pos-10, pos_sensor + 2)
    ejes.set_aspect('equal', 'box')
    ejes.legend()
    
    ang = np.linspace(0, 2*np.pi, 100)
    for H in haces: 
        print(H.rayos[0,0])
        mini = min(H.haz[-1][:,0])
        r = (max(H.haz[-1][:,0])-mini)/2
        if H == haces[0]: 
            ejes2.plot(r*np.cos(ang), (mini+r)*np.ones_like(ang)+r*np.sin(ang), H.color, label='Haz a $\infty$ mm')
        else: ejes2.plot(r*np.cos(ang), (mini+r)*np.ones_like(ang)+r*np.sin(ang), H.color, label=f'Haz a {round(pos-H.posiciones[0])} mm')
    ejes2.set_xlim(-2, 5)
    ejes2.set_title(f'Desenfoque a f{fn}')
    ejes2.legend()
    #plt.show()

'''
Construir  sistemas ópticos más complejos, 
y obtener sus matrices para conseguir puntos y distancias focales
'''
d_final = 73.804
nt = 100
t=np.linspace(0,1,nt)


l = 8*(1-np.exp(-1.17*t))#t**(0.4)#np.linspace(0,12,100)
k = 50.4*(1-np.exp(-4.45*t))#np.linspace(0, 45, nt)
factor = 3
pos=2000
focus1=0#-2.1985
focus = 0
n1 = 1.713
n2 = 1.64831
#H = np.linspace(19, 19, nt) #f2.8
H = np.linspace(2.5, 2.5, nt) #f5.6

Diff = []


rayos_prolongados = False
'''
fig = plt.figure(figsize=(16,8))
ejes = fig.add_subplot(1,2,1, figsize=(10, 5))
ejes2 = fig.add_subplot(1,2,2, figsize=(3, 5))
'''
fig, axs = plt.subplots(1,2,figsize=(16,9), gridspec_kw={'width_ratios': [7, 4]})
ejes = axs[0]
ejes2 = axs[1]


for i in  [90, ]:# [25, 42, 90]:  #np.linspace(25, 99, 74, dtype=int):#range(25, 100, 1): #[25,]:
    
    sist = [#Grupo 1
            Intercara(n0=1,  n1=n1,R01= -88.702, config={'h':28, 'pos':pos + 6.3+focus1}),
            Intercara(n0=n1,  n1=1,R01= np.inf, config={'h':28, 'pos':pos+6.78+focus1}),
            
            Intercara(n0=1,  n1=n1,R01= -37.399, config={'h':22, 'pos':pos +6.78 + 7.50+focus1}),
            Intercara(n0=n1,  n1=1,R01= -61.587, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 +focus1}),
            
            
            Intercara(n0=1,  n1=n2,  R01= -78.525, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 + 0.4+focus}),
            Intercara(n0=n2,  n1=1,R01= -31.585, config={'h':16, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 1+focus}),
            
            Diafragma(pos=2025.5, y= 2.06,#16 para f2.8, 2.06 para f22, 8.1 para f5.6
                      ymax=22),
            
            #Grupo 2
            Intercara(n0=1,  n1=n2,R01= 35.134, config={'h':16, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 7}),
            Intercara(n0=n2,  n1=n1,R01= np.inf, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 8 + 5}),
            Intercara(n0=n1,  n1=1,R01= 46.274, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 8 + 6 }),
            
            Intercara(n0=1,n1=n1,  R01= np.inf, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 8 + 6 + 1 + 4.5}),
            Intercara(n0=n1,  n1=1,R01= 80, config={'h':22, 'pos':pos + 6.78 + 0.19 + 8.72 + 2.71 + 1.94 + 8 + 6 + 1 + 5.5 + 1.62}),
            #Seguir por aquí
            
            ]
    
    pos_sensor = sist[-1].pos + 73.804#71.246 #Montura Sony E

    
    #H1 = Haz([0,28.5], np.linspace(-0.023,-0.005,3), config={'color':'tab:cyan'}) #f2.8
    #H1 = Haz([0,28.5], np.linspace(-0.0185,-0.0095,3), config={'color':'tab:cyan'}) #f5.6
    H1 = Haz([0,28.5], np.linspace(-0.0153,-0.01285,3), config={'color':'tab:cyan'}) #f22
    
    H2 = Haz([0, np.linspace(-H[i], H[i], 3)], 0, config={
        'color': 'y'})
    
    H3 = Haz([0,0], np.linspace(-0.1,0.1,2), config={
        'color':'y'})
    
    #H4 = Haz([300,-133], np.linspace(0.0665,0.088,3), config={'color':'tab:purple'}) #f2.8
    #H4 = Haz([300,-133], np.linspace(0.0715,0.0825,3), config={'color':'tab:purple'}) #f5.6
    H4 = Haz([300,-133], np.linspace(0.07565,0.0784,3), config={'color':'tab:purple'}) #f22
    
    H5 = Haz([50, np.linspace(-85, -80, 5)], 0.17, config={
        'color': 'g'})
    
    haces = [H2, H1, H4]#H1, H4] #Estudiar la distancia de enfoque
    

    ejes.cla()
    ejes2.cla()
    simulacion(sist, haces, ejes, i)
    
    plt.pause(0.05)
    '''
print(max(Diff)-min(Diff))
print(abs(Diff[0]-Diff[-1]))
fig2 = plt.figure(figsize=(16,8))
ejes2 = fig2.add_subplot(1,1,1)
ejes2.plot(np.linspace(25,99,74), Diff)
plt.show()
'''