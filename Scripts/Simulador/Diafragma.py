# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 10:13:07 2022

@author: danie
"""
import numpy as np
import matplotlib.pyplot as plt

class Diafragma:
    def __init__ (self, pos, ymax, y):
        self.pos = pos
        self.ymax = ymax
        self.h = y
        self.M = np.identity(2)
        self.n0=1
        self.n1=1
        self.n01=1
    
    def dibuja(self, ejes):
        
        ejes.plot(self.pos*np.ones(20), np.linspace(self.h, self.ymax, 20), 'k')
        ejes.plot(np.linspace(self.pos-0.5, self.pos+0.5, 20), self.h*np.ones(20), 'k')
        
        ejes.plot(self.pos*np.ones(20), np.linspace(-self.h, -self.ymax, 20), 'k')
        ejes.plot(np.linspace(self.pos-0.5, self.pos+0.5, 20), -self.h*np.ones(20), 'k')