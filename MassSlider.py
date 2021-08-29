

import numpy as np
import gls
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from PyAstronomy.pyTiming import pyPeriod


#SECTION 1: READING DATA
datadir='RV_data/'

data1=datadir+'HD59967.vels'                                                            #File name of first data file
data2=datadir+'HD76653.vels'
data3=datadir+'HD115820.vels'                                                           #Has very few entries
                                                                                      #Has very few entries

star = data1
data = np.genfromtxt(star)                                                     #Reading the file
time = data[:,0]                                                               #Reading time from the file
time = time-time[0]                                                            #Resetting the origin of time to beginning of the observation period
rv = data[:,1]                                                                 #Radial Velocity measured
rv = rv-np.median(rv)                                                          #Subtracting the star's center of mass radial velocity and only keeping the oscillations
error = data[:,2]                                                              #Reading error values off the file

#SECTION 2: CALCULATING PARAMETERS
diff = np.zeros(len(time))                                                     #List for time periods between two successive measurements
for i in range(len(time)-1):
    diff[i]=time[i+1]-time[i]
f_max = 1/np.median(diff)                                                            #Maximum frequency that could be safely extracted
f_min = 1/(max(time)-min(time))                                                #Minimum frequency that could be certainly extracted
res = 500                                                                      #The resolution/least count of frequency
farray = np.linspace(f_min,f_max,res)                                          #The freq range over which we will calculate the periodogram

prim = pyPeriod.Gls((time,rv))
fapp = prim.powerLevel(0.1)

#SECTION 3: DEFINING FUNCTIONS
class pgm:                                                                     #pgm stands for periodogram
    def _init_(self,f,mp,ms,d):                                                #Initialising the class
        self.f=f
        self.mp=mp
        self.ms=ms
        self.d=d
    
    def A(ms,mp,f):
        return 20948.82*f**(1/3)*mp*(mp+ms)**(-2/3)                            # A stands for Amplitude of the RV oscillation. The constant is calculated so that you can input f in 1/day and masses in Mjups to get A in m/s
    
    def FT(f,mp,d):                                                            # Resturns the power values of the Fourier Transform of input signal+ added sine.  f=freq, mp=mass planet, d=phase
        rv1=rv+pgm.A(1047.94,mp,f)*np.sin(2*np.pi*time*f+d)                    #Adding a simulated signal to our RV values
        gls2=pyPeriod.Gls((time,rv1,error),fbeg=f_min,fend=f_max,freq=farray)  #Taking Lomb Scargle Periodogram
        return gls2.power
    
    def fap(f,mp,d,fapv):                                                      # Resturns the power values of the Fourier Transform of input signal+ added sine. f=freq, mp=mass planet, d=phase
        rv1=rv+pgm.A(1047.94,mp,f)*np.sin(2*np.pi*time*f+d)                    #Adding a simulated signal to our RV values
        gls2=pyPeriod.Gls((time,rv1,error),fbeg=f_min,fend=f_max,freq=farray)  #Taking Lomb Scargle Periodogram
        fapvalues=np.zeros(res)
        for i in range(0,res):
            fapvalues[i]=gls2.powerLevel(fapv)                                 #power level gives the power threshold for the given value for FAP
        return fapvalues
    
#SECTION 4:  PLOTTING

"""Below is the code that creates a dynamic plot, you could move the sliders and change parameteres to see an updated graph immediately"""


fig, ax = plt.subplots()                                                       #Standard plt declarations 
plt.subplots_adjust(left=0.13,bottom=0.25)

p,=plt.plot(farray,pgm.FT((f_max+f_min)/2,0,0))                                #farray is the frequency array, the x-axis. pgm.FT 
q,=plt.plot(farray,pgm.fap((f_max+f_min)/2,0,0,0.01))
#r,=plt.plot(farray,pgm.windowf())
plt.xlabel("Frequency spectrum")
plt.ylabel("Power")
plt.title(star)
ax.set_ylim([0,1])
axfreq = plt.axes([0.13, .12, 0.775, 0.03])
axamp = plt.axes([0.13, .07, 0.775, 0.03])
axphase= plt.axes([0.13, .02, 0.775, 0.03])

sldrf = Slider(axfreq, 'Planet Freq', f_min, f_max, valinit=(f_max+f_min)/2,valstep=0.01)
sldrmp = Slider(axamp, 'Planet Mass', 0, 2, valinit=0,valstep=0.01)
sldrd = Slider(axphase, 'Phase',0, 2*np.pi, valinit=0,valstep=0.05)

def update(val):
    f = sldrf.val
    mp = sldrmp.val
    d= sldrd.val
    p.set_ydata(pgm.FT(f,mp,d))
    #q.set_ydata(pgm.fap(f,mp,d,0.1))
    fig.canvas.draw_idle()

sldrf.on_changed(update)
sldrmp.on_changed(update)
sldrd.on_changed(update)

plt.show()
