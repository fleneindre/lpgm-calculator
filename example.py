# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 12:54:02 2021

@author: Francois LN @NoneType1
"""

import scipy.io
import numpy as np
from dvg_ringbuffer import RingBuffer
import matplotlib.pyplot as plt
from LPGMCalculator import LPGMCalculator

mat = scipy.io.loadmat('matsuki-20210213.mat')

accel = mat.get('acc')

# Sample frequency
sampleRate = 100

lpgmCalculator = LPGMCalculator(sampleRate)


accPx = RingBuffer(30*sampleRate, dtype=np.float64)
accPx.extend(np.zeros((30*sampleRate)))
accPy = RingBuffer(30*sampleRate, dtype=np.float64)
accPy.extend(np.zeros((30*sampleRate)))
accPz = RingBuffer(30*sampleRate, dtype=np.float64)
accPz.extend(np.zeros((30*sampleRate)))

i = 0

t = np.linspace(-29.99,0,3000)

while True:

    rawAcc = accel[i,:]
    
    LPGM = lpgmCalculator.update(rawAcc)
    
    accF = lpgmCalculator.getFilteredAcceleration()
    accPx.append(accF[0])
    accPy.append(accF[1])
    accPz.append(accF[2])

   
        
    if i % 100 == 0:
        
        maxSva30 = lpgmCalculator.getMaxSva30()
        vectPx = accPx[:]
        vectPy = accPy[:]
        vectPz = accPz[:]
        PGA = np.max(np.sqrt(vectPx[-100:]**2+vectPy[-100:]**2+vectPz[-100:]**2))
        maxSva = lpgmCalculator.getMaxSva()
        
        valMax = np.max([np.max(np.abs(vectPx)),np.max(np.abs(vectPy)),np.max(np.abs(vectPz))])
        
        strTitle = '福島市松木町   PGA {0:.2f} gal   Sva {1:.2f} cm/s  LPGM {2:d} ({3:.2f} cm/s)'.format(PGA,np.max(maxSva[:100]), LPGM,maxSva30)
        
        plt.clf()
        plt.subplot(3,1,1)
        
        plt.title(strTitle, fontname='Meiryo')
        plt.plot(t,vectPx)
        plt.ylim(-valMax-5, valMax+5)
        plt.margins(x=0)
        plt.ylabel('Acc NS [gal]')
        plt.grid()
        
        plt.subplot(3,1,2)
        plt.plot(t,vectPy)
        plt.ylabel('Acc EW [gal]')
        plt.margins(x=0)
        plt.ylim(-valMax-5, valMax+5)
        plt.grid()
        
        plt.subplot(3,1,3)
        plt.plot(t,vectPz)
        plt.ylabel('Acc UD [gal]')
        plt.xlabel('Time [s]')
        plt.margins(x=0)
        
        plt.ylim(-valMax-5, valMax+5)
        plt.grid()
       
        print(strTitle)
        
     
       
        plt.pause(0.05)
    
   
    i += 1
    
    if i == 30000:
        break
plt.show()
