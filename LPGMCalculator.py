# -*- coding: utf-8 -*-
"""
Created on Sun Dec 19 22:33:50 2021

@author: François Le Neindre / @Nonetype1 

このクラスは、気象庁の長周期地震動階級をリアルタイムの生加速度データから計算するように設計されています。

また、水平成分の絶対速度応答スペクトルを取得することもできます。これは、0.2秒刻みで1.6秒から7.8秒の範囲です。

使用した方法は、[Nigam, Navin C. and Jennings, Paul C. (1969) Calculation of 
response spectra from strong-motion earthquake records.]からのものです。

/！\読み取り値は、加速度計が建物の1階、または地面に接続された最も水平な面に配置されている場合に
のみ有効です。
サンプルレートは、時間の経過とともに一定である必要があります。
加速度入力はフィルタリングされていない必要があります（重力補正なし）


This class is designed to compute the Japan Meteorological Agency Long-Period 
Ground Motion class, from real-time raw acceleration data.
It is also possible to retrieve the Absolute Velocity Response Spectrum of 
horizontal components, spanning between 1.6s and 7.8s with 0.2s increments.

The method used is from [Nigam, Navin C. and Jennings, Paul C. (1969) 
Calculation of response spectra from strong-motion earthquake records.]

/!\ The readings are only valid if the accelerometer is placed on the first 
floor of any building or on a solid surface tied to the ground as level as 
possible. The sample rate has to be constant over time.
The acceleration input has to be unfiltered (ie. without gravity compensation)

"""

import numpy as np
from dvg_ringbuffer import RingBuffer
from scipy import signal

class LPGMCalculator:
    
    def __init__(self,sampleRate):
        
        self.sampleRate = sampleRate
        # Sample time
        
        self.sampleTime = 1/sampleRate

        # High-Pass Filter coefficients (cut-off at 20 seconds)
        baf = signal.butter(2, 0.05, 'hp', fs=sampleRate, output='ba')
        self.bf = baf[0]
        self.af = baf[1]

        # Sva periods array (1.6s to 7.8s with 0.2s increments)
        self.Np = 32
        self.periods = np.linspace(1.6,7.8,self.Np)
        self.Sva = np.zeros((self.Np))
        

        # Damping factor
        beta = 0.05;

        # A and B matrix collection initialization
        self.A = np.zeros((2,2,self.Np))
        self.B = np.zeros((2,2,self.Np))

        # x matrix collection initialization
        self.xi = np.zeros((2,2,self.Np))

        self.init = True

        for j in range(0, self.Np):
            
            w = 2*np.pi/self.periods[j]
            SW = np.sin(w*np.sqrt(1-beta**2)*self.sampleTime)
            CW = np.cos(w*np.sqrt(1-beta**2)*self.sampleTime)
            E = np.exp(-beta*w*self.sampleTime)
            
            self.A[0,0,j] = E*(beta/np.sqrt(1-beta**2)*SW+CW)
            
            self.A[0,1,j] = SW*E/(w*np.sqrt(1-beta**2))
            
            self.A[1,0,j] = -SW*E*w/(np.sqrt(1-beta**2))
            
            self.A[1,1,j] = E*(-beta/np.sqrt(1-beta**2)*SW+CW)
            
            self.B[0,0,j] = E*(((2*beta**2-1)/(w**2*self.sampleTime)+beta/w)*SW/(w*np.sqrt(1-beta**2))+\
                          ((2*beta)/(w**3*self.sampleTime)+1/w**2)*CW)-2*beta/(w**3*self.sampleTime)
            
            self.B[0,1,j] = -E*(((2*beta**2-1)/(w**2*self.sampleTime))*SW/(w*np.sqrt(1-beta**2))+\
                ((2*beta)/(w**3*self.sampleTime))*CW)+2*beta/(w**3*self.sampleTime)-1/w**2
            
            self.B[1,0,j] = E*(((2*beta**2-1)/(w**2*self.sampleTime)+beta/w)*(CW-beta*SW/np.sqrt(1-beta**2))-\
                ((2*beta)/(w**3*self.sampleTime)+1/w**2)*(w*np.sqrt(1-beta**2)*SW+beta*w*CW))+1/(w**2*self.sampleTime)
            
            self.B[1,1,j] = -E*(((2*beta**2-1)/(w**2*self.sampleTime))*(CW-beta*SW/np.sqrt(1-beta**2))-\
                ((2*beta)/(w**3*self.sampleTime))*(w*np.sqrt(1-beta**2)*SW+beta*w*CW))-1/(w**2*self.sampleTime);

        # Buffers initialization
        self.accH = np.array([0,0,0])
        self.accH_1 = np.array([0,0,0])
        self.accH_2 = np.array([0,0,0])

        self.accHF = np.array([0,0,0])
        self.accHF_1 = np.array([0,0,0])

        self.maxSva = RingBuffer(30*sampleRate, dtype=np.float64)
        self.maxSva.extend(np.zeros((30*sampleRate)))

        self.vel = np.array([0,0,0])

        self.acc0 = np.array([0,0,0])  
        
        self.LPGM = 0
        
        
    # Update the LPGM, filtered 3D acceleration, 3D velocity
    # Absolute Velocity Response Spectrum (Sva) of combined horizontal components
    #
    # Input : current raw acceleration components (np.array([ax,ay,az])
    # By convention, the 2 first components have to be horizontal 
    #
    # Output : LPGM value taken from the maximum Sva value of the last 30 seconds
    
    
    def update(self, rawAcceleration):

        if self.init:
            # Store the first value of raw acceleration to compensate the continuous component
            self.acc0 = rawAcceleration
            self.init = False

                
        # Shift filter registers
        self.accH_2 = self.accH_1
        self.accH_1 = self.accH
        self.accH = rawAcceleration - self.acc0
        
        # Filter the acceleration components
        accHFi = self.bf[0]*self.accH + self.bf[1]*self.accH_1 + self.bf[2]*self.accH_2 \
            - self.af[1]*self.accHF - self.af[2]*self.accHF_1
        
        self.accHF_1 = self.accHF
        self.accHF = accHFi
        
        # Compute current velocity
        self.vel = self.vel + (self.accHF_1+self.accHF)*self.sampleTime/2;

        for j in range(0, self.Np):
            
            # Compute current oscillator horizontal relative velocity and displacement
            self.xi[:,:,j] = np.matmul(self.A[:,:,j],self.xi[:,:,j]) + \
                np.matmul(self.B[:,:,j],np.array([self.accHF_1[0:2],self.accHF[0:2]]))
            
            # Store the current oscillator absolute velocity norm
            self.Sva[j] = np.sqrt((self.xi[1,0,j]+self.vel[0])**2+(self.xi[1,1,j]+self.vel[1])**2)
            
        
        # Store the maximum Sva value
        self.maxSva.appendleft(np.max(self.Sva))
        
        # Compute JMA LPGM class from the last 30sec of Sva data 
        self.maxSva30 = np.max(self.maxSva)
        
        if self.maxSva30 < 5:
            self.LPGM = 0
        elif self.maxSva30 < 15:
            self.LPGM = 1
        elif self.maxSva30 < 50:
            self.LPGM = 2
        elif self.maxSva30 < 100:
            self.LPGM = 3
        else:
            self.LPGM = 4
        return self.LPGM

    
    def getSva(self):
        # Returns Absolute Response Velocity Spectrum from 1.6 to 7.8s with 0.2s increment
        return self.Sva
    
    def getFilteredAcceleration(self):
        # Returns the High-Pass filtered acceleration (recommended for PGA computation)
        return self.accHF
    
    def getVelocity(self):
        # Returns the current 3D velocity
        return self.vel
    
    def getMaxSva30(self):
        #Returns the maximum of the maximum Sva values of the last 30 seconds
        return self.maxSva30
    
    def getMaxSva(self):
        #Returns the current maximum Sva value
        return self.maxSva