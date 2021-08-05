#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 16:19:30 2020

@author: yashodhan
"""

from matplotlib.widgets import Slider  # import the Slider widget
import gls # generalised lomb scargle periodogram (Fourier Transform with wavelength as X axis) module
import numpy as np
import matplotlib.pyplot as plt

data1='HD59967.vels' #File name of first data file
data2='HD76653.vels'
data3='HD115820.vels'

data = np.genfromtxt(data2) #Reading the file
time= data[:,0]
time=time-time[0] #Resetting the origin of time to beginning of the observation period
rv=data[:,1] #Radial Velocity measured
rv=rv-np.median(rv) #Subtracting the star's center of mass radial velocity and only keeping the oscillations
error=data[:,2] #Error in RV measurement


diff=[] #List for time periods between two successive measurements
for i in range(len(time)-1):
    diff.append(time[i+1]-time[i])
f_max=1/max(diff) #Maximum frequency of the sine wave hidden in the data
f_low=1/(max(time)-min(time)) #Minimum frequency that could be certainly extracted
farray=np.linspace(f_low,f_max,500)
pgm=gls.Gls((rv,time,error),fbeg=f_low,fend=f_max,freq=farray)
#pgm.plot()
#fo=pgm.freq
po=pgm.power

def plot(rv,time,error):
    #pgm=gls.Gls((rv,time,error),fbeg=f_low,fend=f_max,freq=farray)
    plt.plot(pgm.freq,pgm.power)
    plt.show()
    plt.plot(time,rv,"o")
    pgm.plot()

def mixer(rv,time,error,f,A,d):
    rv=rv+A*np.sin(2*np.pi*time*f+d)
    pgm1=gls.Gls((rv,time,error),fbeg=f_low,fend=f_max,freq=farray)
    power=pgm1.power+po
    plt.plot(pgm1.freq,pgm1.power)
    plt.show()
    plt.plot(pgm.freq,po)
    plt.show()
    plt.plot(pgm1.freq,power)
    
mixer(rv,time,error,0.3,10,0)
#pgm.plot()