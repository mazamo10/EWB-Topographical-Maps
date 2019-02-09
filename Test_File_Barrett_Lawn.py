#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 12:49:50 2018

@author: matthew
"""
import sys
sys.path.insert(0, '/Users/matthew/Documents/EWB/Functions')

import numpy as np
import PolartoCartesian as pol2cart
import contours as con
import matplotlib.pyplot as plt

filename='Barrett Lawn 1.xlsx'
points=pol2cart.PolartoCartesian(filename,sheetnumber=2)

x=points[:,0]
y=points[:,1]
z=points[:,2]
#%% 2D Point Plot

fig1=plt.figure()
ax = fig1.add_subplot(111)
ax.scatter(x,y,c='k')


for X, Y, Z in zip(x, y, np.round(z,2)):
    # Annotate the points 5 _points_ above and to the left of the vertex
    ax.annotate('{}'.format(Z), xy=(X,Y), xytext=(-5, 5), ha='right',
                textcoords='offset points',size=8)


ax.set_xlabel('x axis')
ax.set_ylabel('y axis')
plt.gca().set_aspect('equal', adjustable='box')
fig1.savefig("Barrett_Lawn2_pts.pdf", bbox_inches='tight')

plt.show()

#%%
x=np.reshape(x,(np.size(x),1))
y=np.reshape(y,(np.size(y),1))
z=np.reshape(z,(np.size(z),1))
z=np.round(z)
con.Contours(x,y,z,1/2,"Contour_map_BL2.pdf")
