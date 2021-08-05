#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 05:19:54 2020

This program calculates 2D detection maps with gradient 
More details will be added later, please email me if clarifications are needed.

@author: yashodhan
"""

import numpy as np
import gls
import matplotlib.pyplot as plt
#from matplotlib.widgets import Slider
from PyAstronomy.pyTiming import pyPeriod
#from mpl_toolkits.mplot3d import Axes3D
import time
import os
from datetime import datetime

#SECTION 1: READING DATA
data1='HD59967.vels'                                                            #File name of first data file
data2='HD76653.vels'
data3='HD115820.vels'                                                           #Has very few entries

star=data1
data = np.genfromtxt(star)                                                      #Reading the file
times= data[:,0]                                                                #Reading time from the file
times=times-times[0]                                                            #Resetting the origin of time to beginning of the observation period
rv=data[:,1]                                                                    #Radial Velocity measured
rv=rv-np.median(rv)                                                             #Subtracting the star's center of mass radial velocity and only keeping the oscillations
error=data[:,2]                                                                 #Reading error values off the file

#SECTION 2: CALCULATING PARAMETERS
diff=np.zeros(len(times))                                                       #List for times periods between two successive measurements
for i in range(len(times)-1):
    diff[i]=times[i+1]-times[i]
    
f_max=1/max(diff)                                                               #Maximum frequency that could be safely extracted
f_min=1/(max(times)-min(times))                                                 #Minimum frequency that could be certainly extracted
res = 500                                                                       #The resolution/least count of frequency
farray = np.linspace(f_min,f_max,res)
farray=np.linspace(f_min,f_max,res)                                             #The freq range over which we will calculate the periodogram

prim=pyPeriod.Gls((times,rv))                                                   #This is the original periodogram 
fapp=prim.powerLevel(0.01)                                                      # this is the false alarm probability pertaining to an error of 1%
ms=1047.94                                                                      # Mass of the star, currently set to 1 Msun in Mjup
#SECTION 3: DEFINING FUNCTIONS
class pgm:                                                                      #pgm for periodogram
    def _init_(self,f,mp,ms,d):
        self.f=f
        self.mp=mp
        self.ms=ms
        self.d=d
    
    def A(ms,mp,f):
        return 20948.82*f**(1/3)*mp*(mp+ms)**(-2/3)                             #The constant is calculated so that you can input f in 1/day and masses in Mjups
    
    def FT(f,mp,d):                                                             # Resturns the power values of the Fourier Transform of input signal+ added sine
        rv1=rv+pgm.A(ms,mp,f)*np.sin(2*np.pi*times*f+d)                     #Adding a simulated signal to our RV values, 1047 is solar mass in terms of jupiter masses
        gls2=pyPeriod.Gls((times,rv1,error),fbeg=f_min,fend=f_max,freq=farray)   #Taking Lomb Scargle Periodogram
        return gls2.power
    
    def fap(f,mp,d,fapv):                                                       # Resturns the power values of the Fourier Transform of input signal+ added sine
        rv1=rv+pgm.A(ms,mp,f)*np.sin(2*np.pi*times*f+d)                     #Adding a simulated signal to our RV values
        gls2=pyPeriod.Gls((times,rv1,error),fbeg=f_min,fend=f_max,freq=farray)   #Taking Lomb Scargle Periodogram
        fapvalues=np.zeros(res)
        for i in range(0,res):
            fapvalues[i]=gls2.powerLevel(fapv)                                  #power level gives the power threshold for the given value for FAP
        return fapvalues
    
    def fappowerlevel(f,mp,d,fapv):                                             # Resturns the power values of the Fourier Transform of input signal+ added sine
        rv1=rv#+pgm.A(1047.94,mp,f)*np.sin(2*np.pi*times*f+d) #Adding a simulated signal to our RV values
        gls2=pyPeriod.Gls((times,rv1),fbeg=f_min,fend=f_max,freq=farray)
        return gls2.powerLevel(fapv)
    
    def FTu(f,A):                                                               # Power averaged over phase THIS FUNCTION IS A RELIC, NOT UPDATED DON'T USE 
        power=np.zeros(res)
        for i in range(100):
            ph=i*2*np.pi/100
            rv1=rv+A*np.sin(2*np.pi*times*f+ph)
            gls3=gls.Gls((times,rv1,error),fbeg=f_min,fend=f_max,freq=farray)
            power=power+gls3.power
        return power/100

    def windowf():
        ones=np.zeros(len(times))
        for i in range(len(times)):
            ones[i]=2.0
        window=pyPeriod.Gls((times,np.zeros(len(times))),fbeg=f_min,fend=f_max,freq=farray)
        return window.power



#SECTION 4:  SLIDER PLOT (It is made into a big comment, don't use it together, it might hang your computer)
"""
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
"""
#SECTION 5: DETECTION LIMITS MAPS
s1=time.time()
zres=36                                         # Resolution in the third axis, so the resolution of the probability. zres=10 would mean detection probability can take values from 0 to 1 in steps of 0.1
xyres=100                                       # Resolution/number of divisions of the x and y axis, of frequency
bandwidth=(f_max-f_min)/100                     # This is the half-bandwidth around the input frequency, too wide and you get false detections in the maps 
mplower=0                                       # Lowest mass for the fake signal
mpupper=1                                       # Upper limit for the mass for the fake signal
fakemass=np.linspace(mplower,mpupper,xyres)     # This is the mass space
fakefreq=np.linspace(f_min,f_max,xyres)         # This is the frequency space
fakemass, fakefreq=np.meshgrid(fakemass, fakefreq) #Making them both 2D
rows, cols = (xyres, xyres) 
prob = [[1 for i in range(cols)] for j in range(rows)] # Detection Probability Matrix
peakpow=[[0 for i in range(cols)] for j in range(rows)] # Peak power values
amp=[[0 for i in range(cols)] for j in range(rows)]     #Amplitude values

for i in range(xyres):
    for j in range(xyres):
        amp[i][j]=pgm.A(ms,fakemass[j][j],fakefreq[i][i]) #amplitude values calculated here. They are not used in calculation, only to be recorded in data.
        
def detprob(mp,f):                               # Indices are inputs Calculates the fraction of phase for which the periodogram signal corresponding to the added signal exceeds the False Alarm Probability
    fband=[]                                    # fband is the approximate frequncy band around the input frequency
                 
    powvals=[]                                  # These will be the respective power values for the chosen band
    dp=0                                        # dp stands for detection probability 
    for i in range(0,res):
        if farray[i]>f-bandwidth and farray[i]<f+bandwidth:
            fband.append(i)                     # i is appended and not farray[i] because the power array follows the same indexing.
     
    for i in range(zres):
        powar=pgm.FT(f,mp,i*2*np.pi/zres)
        for j in fband:
            powvals.append(powar[j])
        if max(powvals)>fapp:                   # if the peak of the powervalue exceeds the false alarm probability level then it is a detection
            dp+=1
    return dp/zres, max(powvals)

def execute():                                                                  #This function is created so you can avoid the main calculation when testing other parts of the code             
    """This function uses while condition and for loop to cover
    the mass-frequency field. While condition stops the calculations when 
    we have probability=1 for an entire range of frequency, because
    here on it's the same for all masses all frequencies """
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
    f=open('DATA'+timestamp[-6:-1]+'.txt','a')
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
plt.title("Detection probability map for "+star)
plt.imshow(np.transpose(prob), extent=extent)   # I don't know why I need to take transpose for correct results, see comment on line 173
plt.xlabel("Planet Mass in Jupiter masses")
plt.ylabel("Planet Frequency in 1/days")
plt.colorbar()
plt.show()


writedata(prob, peakpow, fakemass, fakefreq, amp)
    
    
    
    
    
    
    
    
    
    
    
    
    
    



