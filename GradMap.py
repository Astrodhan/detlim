#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 05:19:54 2020

This program calculates 2D detection maps with gradient 
More details will be added later, please email me if clarifications are needed.

@author: yashodhan
"""

import numpy as np
import matplotlib.pyplot as plt
#from matplotlib.widgets import Slider
from PyAstronomy.pyTiming import pyPeriod
#from mpl_toolkits.mplot3d import Axes3D
import time
import os
from datetime import datetime

#SECTION 1: READING DATA

datadir='RV_data/'

data1=datadir+'HD59967.vels'                                                    #File name of first data file
data2=datadir+'HD76653.vels'
data3=datadir+'HD115820.vels'                                                   #Has very few entries
                                                                                #For now, you have to manually add stars like this.

star=data1                                                                      #Change here to change the star
data = np.genfromtxt(star)                                                      #Reading the file
times= data[:,0]                                                                #Reading time from the file
times=times-times[0]                                                            #Resetting the origin of time to beginning of the observation period
rv=data[:,1]                                                                    #Radial Velocity measured
rv=rv-np.median(rv)                                                             #Subtracting the star's center of mass radial velocity and only keeping the oscillations
error=data[:,2]                                                                 #Reading error values off the file

#SECTION 2: CALCULATING PARAMETERS
diff=np.zeros(len(times))                                                       #diff is the array which contains [t(n+1)-t(n)] entries
for i in range(len(times)-1):
    diff[i]=times[i+1]-times[i]
    
f_max=1/max(diff)                                                               #Maximum frequency that could be safely extracted, I have chosen max(diff) instead of min(diff), it should not make much difference in regularly spaced data
f_min=1/(max(times)-min(times))                                                 #Minimum frequency that could be certainly extracted, max(times)-min(times) is the total time-window of observations.
res = 100                                                                       #The resolution/least count of frequency
farray=np.linspace(f_min,f_max,res)                                             #The freq range over which we will calculate the periodogram
falsealarmprob=0.01
prim=pyPeriod.Gls((times,rv,error),fbeg=f_min,fend=f_max,freq=farray)                                                   #This is the original periodogram, prim for primary.
fapp=prim.powerLevel(falsealarmprob)                                                      # this is the false alarm probability pertaining to an error of 1%; if we have a peak in the periodogram above this level it means that there is only 1% chance of that happening randomly.

ms=1047.94                                                                      # Mass of the star, currently set to 1 Msun in Mjup



#SECTION 3: DEFINING FUNCTIONS
class pgm:                                                                      #pgm for periodogram
    def _init_(self,f,mp,ms,d):
        self.f=f
        self.mp=mp
        self.ms=ms
        self.d=d
    
    def A(ms,mp,f):
        return 20948.82*f**(1/3)*mp*(mp+ms)**(-2/3)                             #This gives us the amplitude of RV oscillations in m/s (yashodhan please check) for the combination of particular mass of the planet, sun and frequency. This is calculated through Kepler's third law (actually through newtonian dynamics). The constant 20948.82 is calculated so that you can input f in 1/day and masses in Mjups
    
    def FT(f,mp,d):                                                             # Returns the power values of the Fourier Transform of input signal+ added sine. Input: frequency of the simulation, mass of the planet and phase
        rv1=rv+pgm.A(ms,mp,f)*np.sin(2*np.pi*times*f+d)                         # Adding a simulated signal to our RV values, 1047 is solar mass in terms of jupiter masses
        gls2=pyPeriod.Gls((times,rv1,error),fbeg=f_min,fend=f_max,freq=farray)  #Taking Lomb Scargle Periodogram
        return gls2.power


def emfunc(f):
    """This function returns 1 for normal frequencipppes and 1.616 for those 
    frequencies for which there are already present peaks 
    in the original periodogram (prim)"""     
    i=farray.tolist().index(f)                                                             #It is called the M function, I just named it so because it looks like an 'M' when there are two peaks in the periodogram
    if prim.power[i]<fapp:
        return 1
    else: 
        return 1.5

#SECTION 5: DETECTION LIMITS MAPS
s1=time.time()                                  # Start time for run-time counting
zres=10                                         # Resolution in the third axis, so the resolution of the probability. zres=10 would mean detection probability can take values from 0 to 1 in steps of 0.1
xyres=res                                      # Resolution/number of divisions of the x and y axis, of frequency
bandwidth=(f_max-f_min)/300                     # This is the half-bandwidth around the input frequency, too wide and you get false detections in the maps 
mplower=0                                       # Lowest mass for the fake signal in jupiter masses
mpupper=1                                       # Upper limit for the mass for the fake signal
fakemass=np.linspace(mplower,mpupper,xyres)     # This is the 1D mass space
fakefreq=farray         # This is the 1D frequency space
fakemass, fakefreq=np.meshgrid(fakemass, fakefreq) #Making them both 2D, I don't fully understand how numpy.meshgrid works.
rows, cols = (xyres, xyres) 
prob = [[1 for i in range(cols)] for j in range(rows)] # Detection Probability Matrix, starting all values with 1.
peakpow=[[0 for i in range(cols)] for j in range(rows)] # Peak power values, starting it as a null matrix
amp=[[0 for i in range(cols)] for j in range(rows)]     #Amplitude values

for i in range(xyres):
    for j in range(xyres):
        amp[i][j]=pgm.A(ms,fakemass[j][j],fakefreq[i][i]) #amplitude values calculated here. They are not used in calculation, only to be recorded in data.
        
def detprob(mp,f):                              # Indices are inputs. Calculates the fraction of phase for which the periodogram signal corresponding to the added signal exceeds the False Alarm Probability
    #fband=[]                                    # fband is the approximate frequncy band around the input frequency
                 
    #powvals=[]                                  # These will be the respective power values for the chosen band, a subset of farray.
    dp=0                                        # dp stands for detection probability 
    
    for i in range(zres):
        powar=pgm.FT(f,mp,i*2*np.pi/zres)#-prim.power
        if powar[farray.tolist().index(f)]/emfunc(f)>fapp:   # Everything inside [] is the index of the frequency in the frequency array 
            dp+=1
    return dp/zres, powar[farray.tolist().index(f)]       #This gives you 1 single detection probability. For example '0.5'.

def execute():                                                                  #This function is created so you can avoid the main calculation when testing other parts of the code             
    """This function uses while condition and for loop to cover
    the mass-frequency field. While condition stops the calculations when 
    we have probability=1 for an entire range of frequency, because from
    here on it's the same for all masses and all frequencies """
    i=0                                                                 
    while i<xyres:
        for j in range(xyres):
            prob[i][j]=detprob(fakemass[i][i],fakefreq[j][j])[0]                # the i and j stuff here is really difficult to process, mentally
            peakpow[i][j]=detprob(fakemass[i][i],fakefreq[j][j])[1]
            
        if prob[i].count(1)==len(prob[i]):
            break
        i+=1
            
#SAVING
def writedata(Z,P,X,Y,amp):
    """ This function writes down all the data"""

    timestamp=str(datetime.today())
    os.makedirs('RYPData', exist_ok=True)
    os.chdir('RYPData')
    os.mkdir(timestamp)
    os.chdir(timestamp)
    plt.savefig('map.png')
    f=open('DATA'+timestamp[-6:-1]+'.txt','a')                                  # This is just to create unique names, to avoid files being over-written.
    f.write('Probability Power PlanetMass Frequency Period RVamplitude')
    f.write('\n')
    roff=2       #round off number of decimals
    for i in range(xyres):
        for j in range(xyres):
            f.write(str(round(Z[i][j],roff))+" "+str(round(P[i][j],roff))+" "+str(round(X[i][i],roff))+" "+str(round(Y[j][j],roff))+" "+str(round(1/Y[j][j],roff))+" "+str(round(amp[i][j],roff)))
            f.write('\n')
    f.close()

execute()                                                                       # Hash this if you wish to avoid the detection probability calculation
s2=time.time()                                              
print("Time taken in seconds: ")
print(s2-s1)

#PLOTTING
fig,ax = plt.subplots()
#ax2= ax.twinx()
plt.subplots_adjust(left=0.1, right=0.9)    #Adjusting frame size
extent= mplower,mpupper,f_min,f_max             # Axes for plotting through imshow
plt.title("Detection probability map for "+star+"\nFor False Alarm Probability of 1%")
plt.imshow(np.transpose(prob), extent=extent)   # I don't know why I need to take transpose for correct results, see comment on line 173

plt.xlabel("Mass of the hypothetical planet (Mjup)")
plt.ylabel("Frequency of the hypothetical (1/days)")
plt.colorbar()
plt.show()


writedata(prob, peakpow, fakemass, fakefreq, amp)
    
    
    
    
    
    
    
    
    
    
    
    
    
    



