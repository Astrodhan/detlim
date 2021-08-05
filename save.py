#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 21 12:53:20 2021

This program records data into a datafile.

@author: yashodhan
"""
from datetime import datetime
import os

def writedata():
    path=os.getcwd()
    timestamp=str(datetime.today())
    os.makedirs('RYPData', exist_ok=True)
    os.chdir('RYPData')
    os.mkdir(timestamp)
    os.chdir(timestamp)
    f=open('DATA'+timestamp[-6:-1]+'.txt','a')
    f.write('Probability Power PlanetMass Frequency Period RVamplitude')
    f.write(prob)
    #conn = sqlite3.connect('RYP table '+timestamp+'.db')

"""
def create_table():
    c = conn.cursor()
    c.execute('CREATE TABLE IF NOT EXISTS RecordONE (probability REAL, power REAL, planetmass REAL, frequency REAL, period REAL, rvamplitude REAL )')
    c.close()

def data_entry(prob,power,mp,freq,period,amp):
    c = conn.cursor()
    c.execute("INSERT INTO RecordONE (probability, power, planetmass, frequency, period, rvamplitude) VALUES(1,2,3,4,7,8)")
    conn.commit()
    c.close()
    
    #f.write()
"""
print('lol')
#create_table()
#data_entry()
  
