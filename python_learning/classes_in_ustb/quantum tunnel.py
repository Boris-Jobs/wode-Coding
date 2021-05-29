# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 09:26:28 2020

@author: 13225
"""


import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.constants as phy
from tqdm import tqdm
from scipy import linalg
from matplotlib.animation import FuncAnimation

class QuantumTunnel:
    def __init__(self, potential_wall,
                 mass = 1, hbar=1,
                 xmin=0, xmax=100, ninterval=1600):
        # 将坐标离散化
        self.x = np.linspace(xmin, xmax, ninterval)    
        self.hbar = hbar
        self.mass = mass    
        self.U = np.diag(potential_wall(self.x), 0)
        
        self.wave, self.avgE = self.wave_packet(self.x)
        
        self.Lap = self.laplacian(ninterval)        
        self.H = - hbar**2 / (2*mass) * self.Lap + self.U       
        self.history = {}
        
    def laplacian(self, N):
        '''构造二阶微分算子：Laplacian'''
        dx = self.x[1] - self.x[0]
        return (-2 * np.diag(np.ones((N), np.float32), 0)
            + np.diag(np.ones((N-1), np.float32), 1)
            + np.diag(np.ones((N-1), np.float32), -1))/(dx**2)
    
    def rho(self, psi):
        '''从归一化的波函数计算概率密度'''
        return (np.conjugate(psi) * psi).real
    
    def evolve(self, tfinal=30.0, nt=400):
        t = np.linspace(0, tfinal, nt)
        dt = t[1] - t[0]
        Ut = linalg.expm(-1j * self.H * dt / self.hbar)
        #print('Ut=', Ut)
        psi_list = []
        rho_list = []
        
        psi = np.copy(self.wave)
        psi_list.append(psi)
        rho_list.append(self.rho(psi))
        
        for i in range(nt):
            psi = np.dot(Ut, psi)
            psi_list.append(psi)
            rho_list.append(self.rho(psi))
            
        return t, self.x, psi_list, rho_list
    
    
    def reflect_probability(self, rho_):
        N = len(self.x)
        dx = self.x[1] - self.x[0]
        return np.sum(rho_[:N//2]) * dx 
                
    
    def wave_packet(self, x, kmu=2, ka=20):
        '''kmu: peak momentum
           ka: momentum width parameter
           return the Fourier transformation of 
                  exp(-ka * (k - kmu)^2) * exp(-6j k^2)
        '''
        L = x[-1] - x[0]
        dk = 2 * np.pi / L
        N = len(x)
        k = np.linspace(0, N*dk, N)
        # 动量空间下的高斯波包
        psi_k = np.exp(-ka*(k - kmu)**2) * np.exp(-6j * k**2)
        # 动能期望值
        temp = np.dot(np.diag(k*k, 0)/(2*self.mass), psi_k)
        avgE = np.dot(np.conjugate(psi_k), temp) * dk
        avgE = avgE / self.norm(psi_k, dk)**2
        print('<E>', avgE)
        # 傅里叶变换到坐标空间
        psi = np.fft.ifft(psi_k)
        dx = self.x[1] - self.x[0]
        psi = psi / self.norm(psi, dx)
        return psi, avgE
    
    def norm(self, psi, mesh_size):
        # 归一离散化的波函数
        norm = np.sqrt(np.dot(np.conjugate(psi), psi) * mesh_size)
        return norm   
    
    def plot_wave_packet(self, show_density=True):
        with plt.style.context(['Solarize_Light2']):
            plt.plot(self.x, self.wave.real, label=r'$\psi(x)$')
            if show_density: 
                density = (np.conjugate(self.wave) * self.wave).real
                plt.plot(self.x, density, label='$\psi^*(x)\psi(x)$')
            plt.xlabel(r'$x$')
            plt.legend(loc='best', title="wave packet")
            
    def plot_potential(self):
        with plt.style.context(['Solarize_Light2']):
            plt.plot(self.x, np.diag(self.U))
            plt.ylabel(r'potential')
            plt.xlabel(r'$x$')
            
            
def barrier(x, avgE=2.06, shape="square"):
    '''shape: {square, heavyside, well}'''
    L = x[-1] - x[0]
    if shape == 'square':
        pot = (np.heaviside(x - 0.45 * L, 0.5)-np.heaviside(x - 0.55 * L, 0.5)) * avgE
    elif shape == 'heavyside':
        pot = np.heaviside(x - 0.5 * L, 0.5) * avgE
    elif shape == 'well':
        pot = (np.heaviside(x - 0.55 * L, 0.5)-np.heaviside(x - 0.45 * L, 0.5)) * avgE
    return pot


# matplotlib inline
pot = lambda x: barrier(x, shape='heavyside')
qt = QuantumTunnel(potential_wall = pot)

#matplotlib inline
pot = lambda x: barrier(x, shape='square')
qt = QuantumTunnel(potential_wall = pot)

qt.plot_potential()# 图1

#qt.plot_wave_packet()# 图2

t, x, psi_list, rho_list = qt.evolve()

#matplotlib notebook

def update(i):
    line.set_data(qt.x, rho_list[i])
    text.set_text(r'$t=%.2f$'%t[i])
    return line, text,
    

potential = pot(qt.x)

fig1, ax1 = plt.subplots()
plt.plot(qt.x, potential * 0.08)
line, = plt.plot(qt.x, rho_list[0])
text = plt.text(0, 0.05, '')
plt.grid(ls="--")

plt.ylabel('probability density')

plt.xlabel(r'$x$')
anim1 = FuncAnimation(fig1, update, frames=400, interval=100, blit=True)
anim1.save('C:/Users/13225/Desktop/quantum_tuneling.mp4')
plt.show()   