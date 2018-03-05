#!/usr/bin/env python
'''
Reads data from zadatakData.txt
Data is assumed to have ' ' as delimiter, and be in format:
time flow preassure position velocity
'''

from numpy import loadtxt
import matplotlib.pyplot as plt

inputFile = open('zadatakData.txt')
lines = loadtxt("zadatakData.txt", delimiter=" ", unpack=False)

t = lines[:,0]
m_dot = lines[:,1]
p = lines[:,2]
x = lines[:,3]
v = lines[:,4]

plt.figure('Zadatak')
ax1 = plt.subplot(311)
plt.plot( t, m_dot ,'g' )
plt.title( 'Protok' )
ax2 = plt.subplot(312, sharex=ax1)
plt.plot( t ,p ,'g' )
plt.title( 'Tlak fluida' )
ax3 = plt.subplot(313, sharex=ax1, sharey=ax1)
plt.plot( t ,x ,'g' )
plt.title( 'Pomak klipa' )
plt.show()

print("Drawing graphs.")
