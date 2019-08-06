import matplotlib.pyplot as plt
import numpy as np
from scipy import special

C_OB = 1 #mol/m3
C_RB = 1 #mol/m3
D_O = 1E-9 #m2/s
D_R = 1E-9 #m2/s
F = 96485 #C/mol
R = 8.314 #J / mol. K
pi = 3.14 
A = 1E-4 #1m2
T = 298 #K
epsilon = np.sqrt(D_O / D_R)
color=['black','red','blue','purple','c','grey','green']
C_Rt = []
C_Ot = []
C_R = []
C_O = []
t_i = []
i_a = []
i_c = []
i_ai = []
i_ci = []


def CurrentVsTime(E_ai,E_ci,E_0,n): 
#applied anodic voltage, appllied cathodic voltage, equilibrium voltage, 
#number of electrons
  E_a = []
  E_c = []
  t = np.linspace(0.001, 1, 100)
  E_a.append(E_ai)
  E_c.append(E_ci)
  while np.abs(E_ai-E_0) < 0.1:
    theta_a = np.exp(n * F * (E_ai - E_0) / R / T)
    theta_c = np.exp(n * F * (E_ci - E_0) / R / T)
    i_ai = C_RB * n * F *A * np.sqrt(D_R / pi / t) *theta_a * epsilon / \
    (1 + theta_a * epsilon)
    i_ci = - C_OB * n * F *A * np.sqrt(D_O / pi / t) / \
    (1 + theta_c * epsilon)
    i_a.append(i_ai)
    i_c.append(i_ci)
    E_ai = E_ai + 0.02
    E_ci = E_ci - 0.02
    E_a.append(E_ai)
    E_c.append(E_ci)

  plt.figure()
  
  for i in range (0,len(E_a)-1):
    overpotential = E_a[i] - E_0
    plt.plot(t,i_a[i],'s-',color = color[i],\
      label='Anodic current at overpotential of %.3fV' \
      % overpotential,markevery=50)
  for i in range (0,len(E_c)-1):
    overpotential = E_c[i] - E_0
    plt.plot(t,i_c[i],'*-',color = color[i],\
      label='Cathodic current at overpotential of %.3fV' \
      % overpotential,markevery=50)
  del E_c[:]
  del E_a[:]
  ax = plt.gca()
  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')
  ax.xaxis.set_ticks_position('bottom')
  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')
  ax.xaxis.set_ticks_position('bottom')
  ax.yaxis.set_ticks_position('left')
  ax.spines['bottom'].set_position(('data',0))
  ax.spines['left'].set_position(('data',0))
  plt.ylabel(r'Current_A')
  plt.xlabel(r'time_s')
  plt.legend(loc='best',edgecolor='none')
  plt.savefig('current_vs_time.png',dpi=300)
  plt.show()
  
def ConcentrationProfileReduction(E_c,E_0,n):
#appllied cathodic voltage, equilibrium voltage, number of electrons
  x = np.linspace(0,2E-4,1000)
  theta_c = np.exp(n * F * (E_c - E_0) / R / T)
  t = 0.001
  t_i.append(t)
  while t < 11:
    C_Ot = C_OB - C_OB * special.erfc(x / 2 / np.sqrt(D_O * t)) / \
    (1 + epsilon * theta_c)
    C_Rt = C_OB * epsilon * special.erfc(x / 2 / np.sqrt(D_O * t)) / \
    (1 + epsilon * theta_c)
    C_R.append(C_Rt)
    C_O.append(C_Ot)
    t = t*10
    t_i.append(t)
  plt.figure()
  for i in range(0,len(t_i)-1):
    plt.plot(x,C_O[i],'o-',color = color[i],\
      label='t=%.3fs,C_O(x,t)' % t_i[i],markevery=50)
  for i in range(0,len(t_i)-1):
    plt.plot(x,C_R[i],'*-',color = color[i],\
      label='t=%.3fs,C_R(x,t)' % t_i[i],markevery=50)
  del C_R[:]
  del C_O[:]
  del t_i[:]
  ax = plt.gca()
  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')
  ax.xaxis.set_ticks_position('bottom')
  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')
  ax.xaxis.set_ticks_position('bottom')
  ax.yaxis.set_ticks_position('left')
  ax.spines['bottom'].set_position(('data',0))
  ax.spines['left'].set_position(('data',0))
  plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
  plt.ylabel(r'Concentraion_$mol/m^3$')
  plt.xlabel(r'x_m')
  plt.legend(loc='best',edgecolor='none')
  plt.title("Cathodic voltage:%.2fV; Equilibrium potential:%.2fV" % (E_c,E_0))
  plt.savefig('Concentration_profile_reduction.png',dpi=300)

def ConcentrationProfileOxidation(E_a,E_0,n):
#applied anodic voltage, equilibrium voltage, number of electrons
  x = np.linspace(0,2E-4,1000)
  theta_a = np.exp(n * F * (E_a - E_0) / R / T)
  t = 0.001
  t_i.append(t)
  while t < 11:
    C_Rt = C_RB - C_RB * theta_a * epsilon * special.erfc(x / 2 / \
      np.sqrt(D_R * t)) / (1 + epsilon * theta_a)
    C_Ot = C_RB * theta_a * special.erfc(x / 2 / \
      np.sqrt(D_R * t)) / (1 + epsilon * theta_a)
    C_R.append(C_Rt)
    C_O.append(C_Ot)
    t = t*10
    t_i.append(t)
  plt.figure()
  for i in range(0,len(t_i)-1):
    plt.plot(x,C_O[i],'o-',color = color[i],\
      label='t=%.3fs,C_O(x,t)' % t_i[i],markevery=50)
  for i in range(0,len(t_i)-1):
    plt.plot(x,C_R[i],'*-',color = color[i],\
      label='t=%.3fs,C_R(x,t)' % t_i[i],markevery=50)
  del C_R[:]
  del C_O[:]
  del t_i[:]
  ax = plt.gca()
  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')
  ax.xaxis.set_ticks_position('bottom')
  ax.spines['right'].set_color('none')
  ax.spines['top'].set_color('none')
  ax.xaxis.set_ticks_position('bottom')
  ax.yaxis.set_ticks_position('left')
  ax.spines['bottom'].set_position(('data',0))
  ax.spines['left'].set_position(('data',0))
  plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
  plt.ylabel(r'Concentraion_$mol/m^3$')
  plt.xlabel(r'x_m')
  plt.legend(loc='best',edgecolor='none')
  plt.title("Anodic voltage:%.3fV; Equilibrium potential:%.3fV" % (E_a,E_0))
  plt.savefig('Concentration_profile_oxidation.png',dpi=300)
  plt.show()

    
CurrentVsTime(1.02,0.98,1,1)
ConcentrationProfileReduction(0.95,1,1)
ConcentrationProfileOxidation(1.05,1,1)