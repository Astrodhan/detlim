#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 05:19:54 2020

@author: yashodhan
"""

import numpy as np
import gls
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from PyAstronomy.pyTiming import pyPeriod
from mpl_toolkits.mplot3d import Axes3D

#SECTION 1: READING DATA
data1='HD59967.vels' #File name of first data file
data2='HD76653.vels'
data3='HD115820.vels' #Has very few entries

star=data1
data = np.genfromtxt(star) #Reading the file
time= data[:,0] #Reading time from the file
time=time-time[0] #Resetting the origin of time to beginning of the observation period
rv=data[:,1] #Radial Velocity measured
rv=rv-np.median(rv) #Subtracting the star's center of mass radial velocity and only keeping the oscillations
error=data[:,2] #Reading error values off the file

#SECTION 2: CALCULATING PARAMETERS
diff=np.zeros(len(time)) #List for time periods between two successive measurements
for i in range(len(time)-1):
    diff[i]=time[i+1]-time[i]
f_max=1/max(diff) #Maximum frequency that could be safely extracted
f_min=1/(max(time)-min(time)) #Minimum frequency that could be certainly extracted
farray=np.linspace(f_min,f_max,500) #The freq range over which we will calculate the periodogram

prim=pyPeriod.Gls((time,rv))
fapp=prim.powerLevel(0.1)

#SECTION 3: DEFINING FUNCTIONS
class pgm:                                              #pgm for periodogram
    def _init_(self,f,mp,ms,d):
        self.f=f
        self.mp=mp
        self.ms=ms
        self.d=d
    
    def A(ms,mp,f):
        return 20948.82*f**(1/3)*mp*(mp+ms)**(-2/3) #The constant is calculated so that you can input f in 1/day and masses in Mjups
    
    def FT(f,mp,d): # Resturns the power values of the Fourier Transform of input signal+ added sine
        rv1=rv+pgm.A(1047.94,mp,f)*np.sin(2*np.pi*time*f+d) #Adding a simulated signal to our RV values
        gls2=pyPeriod.Gls((time,rv1,error),fbeg=f_min,fend=f_max,freq=farray) #Taking Lomb Scargle Periodogram
        return gls2.power
    
    def fap(f,mp,d,fapv): # Resturns the power values of the Fourier Transform of input signal+ added sine
        rv1=rv+pgm.A(1047.94,mp,f)*np.sin(2*np.pi*time*f+d) #Adding a simulated signal to our RV values
        gls2=pyPeriod.Gls((time,rv1,error),fbeg=f_min,fend=f_max,freq=farray) #Taking Lomb Scargle Periodogram
        fapvalues=np.zeros(500)
        for i in range(0,500):
            fapvalues[i]=gls2.powerLevel(fapv) #power level gives the power threshold for the given value for FAP
        return fapvalues
    
    def fappowerlevel(f,mp,d,fapv): # Resturns the power values of the Fourier Transform of input signal+ added sine
        rv1=rv#+pgm.A(1047.94,mp,f)*np.sin(2*np.pi*time*f+d) #Adding a simulated signal to our RV values
        gls2=pyPeriod.Gls((time,rv1),fbeg=f_min,fend=f_max,freq=farray)
        return gls2.powerLevel(fapv)
    
    def FTu(f,A): # Power averaged over phase THIS FUNCTION IS A RELIC, NOT UPDATED DON'T USE 
        power=np.zeros(500)
        for i in range(100):
            ph=i*2*np.pi/100
            rv1=rv+A*np.sin(2*np.pi*time*f+ph)
            gls3=gls.Gls((time,rv1,error),fbeg=f_min,fend=f_max,freq=farray)
            power=power+gls3.power
        return power/100

    def windowf():
        ones=np.zeros(len(time))
        for i in range(len(time)):
            ones[i]=2.0
        window=pyPeriod.Gls((time,np.zeros(len(time))),fbeg=f_min,fend=f_max,freq=farray)
        return window.power



#SECTION 4:  PLOTTING

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.13,bottom=0.25)
p,=plt.plot(farray,pgm.FT((f_max-f_min)/2,0,0))
q,=plt.plot(farray,pgm.fap((f_max-f_min)/2,0,0,0.1))
#r,=plt.plot(farray,pgm.windowf())
plt.xlabel("Frequency spectrum")
plt.ylabel("Power")
plt.title(star)
ax.set_ylim([0,1])
axfreq = plt.axes([0.13, .12, 0.775, 0.03])
axamp = plt.axes([0.13, .07, 0.775, 0.03])
axphase= plt.axes([0.13, .02, 0.775, 0.03])

sldrf = Slider(axfreq, 'Frequency', f_min, f_max, valinit=(f_max-f_min)/2,valstep=0.01)
sldrmp = Slider(axamp, 'Planet Mass', 0, 10, valinit=0,valstep=0.01)
sldrd = Slider(axphase, 'Phase',0, 2*np.pi, valinit=0,valstep=0.05)

def update(val):
    f = sldrf.val
    mp = sldrmp.val
    d= sldrd.val
    p.set_ydata(pgm.FT(f,mp,d))
    q.set_ydata(pgm.fap(f,mp,d,0.1))
    fig.canvas.draw_idle()

sldrf.on_changed(update)
sldrmp.on_changed(update)
sldrd.on_changed(update)

plt.show()

#SECTION 5: DETECTION LIMITS MAPS

def detprob(f,mp):                              #Calculates the fraction of phase for which the periodogram signal corresponding to the added signal exceeds the False Alarm Probability
    fband=[]                                    #fband is the approximate frequncy band around the peak
    powvals=[]                                  #These are the power values in corresponding to that band, containing the maxima
    dp=0
    for i in range(0,500):
        if farray[i]>f-(f_max-f_min)/10:
            if farray[i]<f+(f_max-f_min)/10:
                fband.append(i)
     
    for i in range(100):
        powar=pgm.FT(f,mp,i*2*np.pi/100)
        for j in fband:
            powvals.append(powar[j])
        if max(powar)>fapp:
            dp+=1
    return dp/100

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X=np.linspace(0,2,100)
Y=np.linspace(f_min,f_max,100)
X,Y=np.meshgrid(X,Y)
rows, cols = (100, 100) 
Z = [[0 for i in range(cols)] for j in range(rows)] 

for i in range(100):
    for j in range(100):
        Z[i][j]=detprob(Y[i][i],X[j][j])
ax.plot_surface(X,Y,Z)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



#gls1.plot()