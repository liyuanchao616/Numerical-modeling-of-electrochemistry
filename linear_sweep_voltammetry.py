# -*- coding: utf-8 -*-
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy import special
import math


C_OB = 1. #mol/m3
C_RB = 1. #mol/m3
D_O = 1E-9 #m2/s
D_R = 1E-9 #m2/s
F = 96485. #C/mol
R = 8.314 #J / mol. K
pi = 3.14 
A = 1E-4 #1m2
T = 298. #K
epsilon = np.sqrt(D_O / D_R)
color=['black','red','blue','purple','c','grey','green']
S_t = []
F_t = []

def CurrentVsVoltage(E_i,E_0,n,v):
 #applied anodic voltage, equilibrium voltage, number of electrons，scan rate (V/s)
  current = [] # current,I(t)
  voltage = [] # voltage,V(t)
  t = np.linspace(0, 50, 101)
  delta_t = (t[-1]-t[0])/(len(t)-1) #step of time
  theta = np.exp(n * F * (E_i - E_0) / R / T)
  epsilon = np.sqrt(D_O/D_R)
  S_t = np.exp(n*F*(-v*t)/R/T) 
  F_t =  - n*F*A*C_OB * np.sqrt(pi*D_O) / (1. + theta*epsilon*S_t)
  s =  Pfuncion_Midpoint(1,0) * F_t[1]/(2.*np.sqrt(delta_t)) #current at t=t[1]
  voltage.append(E_i)
  voltage.append(E_i-v*delta_t)
  current.append(0) #initial current=0
  current.append(s)
  
  
  for i in range(2,len(t)):
    k = 1
    sum = 0
    while k <= i-1:
      sum = sum + (Pfuncion_Midpoint(i,k-1)-Pfuncion_Midpoint(i,k)) * current[k]
      k = k +1
    s = (1/Pfuncion_Midpoint(i,i-1)) * (F_t[i]/(2.*np.sqrt(delta_t))-sum)
    E = E_i - v*i*delta_t
    current.append(s)
    voltage.append(E)


  plt.figure()
  plt.plot(voltage,current,'-',color = 'red')
  ax = plt.gca()
  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')
  ax.xaxis.set_ticks_position('bottom')
  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')
  ax.xaxis.set_ticks_position('bottom')
  ax.yaxis.set_ticks_position('left')

  plt.ylabel(r'Current_A')
  plt.xlabel(r'voltage_V')
  plt.legend(loc='best',edgecolor='none')
  plt.savefig('current_vs_potential.png',dpi=300)
  plt.show()
  
def Pfuncion_PartialMidpoint(n,i):
  p = (2./3) * ((n-i)**1.5-(n-i-1)**1.5)
  return p
def Pfuncion_Midpoint(n,i):
  p = np.sqrt(n-i-1./2)
  return p   
CurrentVsVoltage(1.1,1,1,0.02)
