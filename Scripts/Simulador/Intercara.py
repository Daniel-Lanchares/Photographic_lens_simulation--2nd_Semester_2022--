# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 10:09:24 2022

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
        

class Diafragma:
    def __init__ (self, pos, ymax, y):
        self.pos = pos
        self.ymax = ymax
        self.y = y
        self.M = np.identity(2)
        self.n0=1
        self.n1=1
        self.n01=1
    
    def dibuja(self, ejes):
        
        ejes.plot(self.pos*np.ones(20), np.linspace(self.y, self.ymax, 20), 'k')
        ejes.plot(np.linspace(self.pos-0.5, self.pos+0.5, 20), self.y*np.ones(20), 'k')