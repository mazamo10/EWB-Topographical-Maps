#%% Site 1 Topographical Map--Matthew Zamora

#%%
import math
import time
import pandas as pd
import pyautogui as auto
from numpy import *
from scipy.optimize import *
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
#%% Opening autoCAD
[p,q]=auto.position()       # remove in final version

time.sleep(2)               # remove in final version
#%%
# Proccess for opening autoCAD for mac computer. Coordinate corresdponds to Spotlight Search icon
# (Magnifine glass in top right of screen). Please use mouse_position.py to determine location on your device.
# Delay used to ensure that app is fully open before the next step occurs.

auto.moveTo(x=1345, y=9, duration=1)
auto.click()
auto.typewrite('autoCAD\n', interval=0.5)  # useful for entering text, newline is Enter
time.sleep(20)

#%%
# Selects format "acad.dwt" for autoCAD file, and opens file format. As mentioned above, use mouse_position.py 
# to determine coorisponding location on your device.
auto.moveTo(x=723, y=478, duration=1)
auto.click()

auto.moveTo(x=1061, y=740, duration=1)
auto.click()

#%%
# Waits for autoCAD to completly open, then moves to command line and changes the units to feet. As mentioned above, use mouse_position.py 
# to determine coorisponding location on your device.
time.sleep(20)
auto.moveTo(x=346, y=792, duration=1)
auto.click()
auto.typewrite('-DWGUNITS', interval=0.2)  # useful for entering text, newline is Enter
auto.typewrite(['enter','2','enter','2','enter','4','enter','y','enter','y','enter','y','enter','y','enter'], interval=0.2)
#%% Refrence values
a=91                                # Distance in feet between A and instrument
b=86                                # Distance in feet between B and instrument
theta_ref_A=330                     # Reference angle in degrees measured by instrument to reference point A
theta_ref_B=359                     # Reference angle in degrees measured by instrument to reference point B
C=abs(theta_ref_B-theta_ref_A)      # Angle in degrees opposite c
C=math.radians(C)                   # Angle in radians opposite c

#%% Refrence Value Calculations
c=math.sqrt(a**2+b**2-2*a*b*math.cos(C))          # Distance in feet between A and B
A=math.acos((a**2-b**2-c**2)/(-2*b*c))            # Angle in radians opposite a
B=math.acos((b**2-a**2-c**2)/(-2*a*c));           # Angle in radians opposite b

#%% Instrument Location
# The location of the instrument is determined by finding the positive
# intersect of two circles of radius a and b. This intersect represents 
# the (x,y) coordinates of the instrument.
def Circle_intersect(circle):
    x=circle[0]
    y=circle[1]
    
    F=empty((2))
    F[0]=math.sqrt(pow(a,2)-pow(x,2))-y
    F[1]=math.sqrt(pow(b,2)-pow((x+c),2))-y
    return F

circle_guess=array([1,1])
    
circle=fsolve(Circle_intersect,circle_guess)

#%% Refrence values
R=matrix([[1000,2000,3000],[1000-c,2000,3000],[1000+circle[0],2000+circle[1],3000]]) # (x,y,z) locations of reference points A and B and the instrument

#%% Excel Data
# Reads excel file containing four columns. The first column corrisponds to height,
# the second to distance measured from the instrument, the third to angles measured by the instrument,
# the last column is a boolean array of 1,0 signifying locations at which the midpoint of creek (or other significant points) are.

survey=pd.read_excel('Basic_Topo1_desk.xlsx')
survey.values

Height=survey['Height']
Height=Height.values
Height=reshape(Height,(size(Height),1))

Distance=survey['Distance']
Distance=Distance.values
Distance=reshape(Distance,(size(Distance),1))

Angle=survey['Angle']
Angle=Angle.values

Key_Feature=survey['Key_Feature']
Key_Feature=Key_Feature.values
Key_Feature=reshape(Key_Feature,(size(Key_Feature),1))

#%%Theta Conversions
#vA new theta is calculated by adjusting the axis so that the x axis lies
# parallel to the line between the reference points A,B. Theta is then used
# to find Cartesian coordinates corresponding to the locations surveyed.

theta=radians([Angle])                  # Converts from degrees to radians
theta=transpose(theta)
theta_ref_A=radians(theta_ref_A)        # Converts reference angle A from degrees to radians
theta=((theta_ref_A)-B)-(theta)         # Converts data so that angle 0 is parallel to the line y=a where a is an arbitrary constant

def Theta_correct(theta):
    
    Theta = empty((size(theta),1)) 

    for i in range(size(theta)):
        if theta[i,0]>pi:
            Theta[i,0]=theta[i,0]-2*pi
        elif theta[i,0]<-pi:
            Theta[i,0]=theta[i,0]+2*pi
        else:
            Theta[i,0]=theta[i,0]
    return Theta

Theta=Theta_correct(theta)

#%% Coordinates
# Uses angle Theta and the distance found between the instrument and the relevant point 
# to determine the cartisian coordinate of that point.

def Finding_Coordinates(Distance,theta):
    Dx=empty((size(Theta),1))
    Dy=empty((size(Theta),1))
    
    for i in range(size(theta)):
        if (Theta[i,0] >= -pi/2 and Theta[i,0] <= pi/2):
            Dy[i,0]=Distance[i,0]*math.sin(Theta[i,0])
            Dx[i,0]=Distance[i,0]*math.cos(Theta[i,0])
        elif (Theta[i,0] > pi/2):
            Dy[i,0]=Distance[i]*math.sin(pi-Theta[i,0])
            Dx[i,0]=-1*(Distance[i,0]*math.cos(pi-Theta[i,0]))
        else:
            Dy[i,0]=-1*(Distance[i,0]*math.sin(pi+Theta[i,0]))
            Dx[i,0]=-1*(Distance[i,0]*math.cos(pi+Theta[i,0]))
    return Dx, Dy

[Dx,Dy]=Finding_Coordinates(Distance,theta)

X=R[2,0]+Dx        # Calculates the x-location of the surveyed points by adding the x-position of the instrument to the x-displacement off of the instrument 
Y=R[2,0]+Dy        # Calculates the y-location of the surveyed points by adding the y-position of the instrument to the y-displacement off of the instrument
Z=3000-Height      # Calculates the z-position at each (x,y) pair
Z=round_(Z)        # Rounds the height value to the nearest interger

data=hstack((X,Y)) # Creates array with two columns corrisponding to (x,y) tuplets


#%% Garpahing points with autoCAD

# Saves as text file with .scr extension for autoCad to read as script file. Script file begins with command
# "_multiple _point" to prompt autocad to plot the following points
Data=savetxt("Data1.scr",data,delimiter=",",header="_multiple _point",comments="") 

# Opens Layer Property Manager
auto.moveTo(x=346, y=792, duration=1)
auto.click()
auto.typewrite('-LAYER', interval=0.2)  # useful for entering text, newline is Enter
auto.typewrite(['enter'], interval=0.2)
auto.typewrite('New', interval=0.2)
auto.typewrite(['enter'], interval=0.2)
auto.typewrite('Survey Points', interval=0.2)
auto.typewrite(['enter'], interval=0.2)
auto.typewrite(['esc'], interval=0.2)

# Opens Layer Selection options
auto.moveTo(x=1326,y=145,duration=1)
auto.mouseDown()
time.sleep(.5)
auto.mouseUp()

# Selects Layer
auto.moveTo(x=None,y=200,duration=1)
auto.mouseDown()
time.sleep(.5)
auto.mouseUp()

auto.moveTo(x=550,y=425,duration=1)
auto.click()

# Moves to and Clicks Command Line, then inputs Survey points
auto.click(x=346, y=792, clicks=1, button='left')
auto.typewrite('SCR', interval=0.2)
auto.typewrite(['enter'], interval=0.2)
auto.moveTo(x=988,y=139,duration=1)
auto.click()
auto.typewrite('Data1.scr', interval=0.2)
auto.typewrite(['enter'], interval=0.2)
auto.moveTo(x=612,y=210,duration=2)
auto.click()
auto.moveTo(x=1034,y=563,duration=2)
auto.click()

time.sleep(2)
auto.moveTo(x=346, y=792, duration=1)
auto.click()
auto.typewrite(['esc'], interval=0.2)
auto.typewrite(['z','enter','e','enter'], interval=0.2)
time.sleep(1)
#%% 2D Point Plot

fig1=plt.figure()
ax = fig1.add_subplot(111)

for i in range(size(X)):
    if Key_Feature[i,0]==1:
        ax.scatter(X[i],Y[i],c='b')
    else:
        ax.scatter(X[i],Y[i],c='k')


ax.set_xlabel('x axis')
ax.set_ylabel('y axis')
plt.gca().set_aspect('equal', adjustable='box')

plt.show()


#%% Delaunay Triangulation

X=reshape(X,(size(X),1))
Y=reshape(Y,(size(Y),1))
Z=reshape(Z,(size(Z),1))

Points=hstack((X,Y))


#%%

tri=Delaunay(Points)

fig2=plt.figure()
ax = fig2.add_subplot(111)

ax.triplot(Points[:,0],Points[:,1],tri.simplices)

ax.plot(Points[:,0], Points[:,1],'ob')


ax.set_xlabel('x axis')
ax.set_ylabel('y axis')
plt.gca().set_aspect('equal', adjustable='box')

plt.show()

#%% Countours

triangles=tri.simplices
Points=hstack((X,Y,Z))

fig3=plt.figure()
plt.gca().set_aspect('equal',adjustable='box')
ax = fig3.add_subplot(111)

# Opens Layer Property Manager
auto.moveTo(x=346, y=792, duration=1)
auto.click()
auto.typewrite('-LAYER', interval=0.2)  
auto.typewrite(['enter'], interval=0.2)
auto.typewrite('New', interval=0.2)
auto.typewrite(['enter'], interval=0.2)
auto.typewrite('Contour Lines', interval=0.2)
auto.typewrite(['enter'], interval=0.2)
auto.typewrite(['esc'], interval=0.2)

# Opens Layer Selection options
auto.moveTo(x=1326,y=145,duration=1)
auto.mouseDown()
time.sleep(.5)
auto.mouseUp()

# Selects Layer
auto.moveTo(x=None,y=200,duration=1)
auto.mouseDown()
time.sleep(.5)
auto.mouseUp()

auto.moveTo(x=550,y=425,duration=1)
auto.click()


for i in range(size(triangles,0)):
    tri_single_a=int(triangles[i,0])
    tri_single_b=int(triangles[i,1])
    tri_single_c=int(triangles[i,2])
    
    step_size=1/2
#%%  Countour Lines for Triangle Side AB  
    step1=int(abs(Points[tri_single_a,2]-Points[tri_single_b,2])/step_size)
    
    if step1==0:
        step_point_ab=array([nan,nan,nan])
        step_point_ab=reshape(step_point_ab,(1,3))
    else:
        step_point_ab=empty((step1,3))
    if all(isnan(step_point_ab))==False:
        for n in range(step1):
            if n==0:
                step_point_ab[n,0]=((1/(2*step1))*Points[tri_single_b,0])+(1-(1/(2*step1)))*Points[tri_single_a,0] # Creates value closest to Points[tri_single_a,0]
                step_point_ab[n,1]=((1/(2*step1))*Points[tri_single_b,1])+(1-(1/(2*step1)))*Points[tri_single_a,1] # Creates value closest to Points[tri_single_a,1]
             
                step_point_ab_lastx=((1/(2*step1))*Points[tri_single_a,0])+(1-(1/(2*step1)))*Points[tri_single_b,0]# Creates value closest to Points[tri_single_b,0]
                step_point_ab_lasty=((1/(2*step1))*Points[tri_single_a,1])+(1-(1/(2*step1)))*Points[tri_single_b,1]# Creates value closest to Points[tri_single_b,1]
             
            elif n==(step1-1):
                step_point_ab[n,0]=step_point_ab_lastx # Assigns x value of last point
                step_point_ab[n,1]=step_point_ab_lasty # Assigns y value of last point    
            else:
                step_point_ab[n,0]=((n/(step1-1))*step_point_ab_lastx)+(1-(n/(step1-1)))*step_point_ab[0,0] # Determines x values of step points between the fisrst and last step point value iterating from points closest to the first value
                step_point_ab[n,1]=((n/(step1-1))*step_point_ab_lasty)+(1-(n/(step1-1)))*step_point_ab[0,1] # Determines y values of step points between the fisrst and last step point value iterating from points closest to the first value
    
            if Points[tri_single_a,2]< Points[tri_single_b,2]: # If point a is lower than point b
                step_point_ab[n,2]=Points[tri_single_a,2]+(n+1)*step_size # Creates z values that increase between a and b
            elif Points[tri_single_b,2]< Points[tri_single_a,2]: # If point b is lower than point a
                step_point_ab[n,2]=Points[tri_single_b,2]+(n+1)*step_size # Creates z values that decrease between a and b
            else: # If z values of points a and b are equal
                step_point_ab[n,2]=Points[tri_single_b,2] # z value is the same as both points
        if Points[tri_single_b,2]< Points[tri_single_a,2]:
            step_point_ab[:,2]=flip(step_point_ab[:,2],0)
        else:
            step_point_ab[:,2]=step_point_ab[:,2]
    else:
         pass
     
#%%  Countour Lines for Triangle Side BC  
    step2=int(abs(Points[tri_single_b,2]-Points[tri_single_c,2])/step_size)
    
    if step2==0:
        step_point_bc=array([nan,nan,nan])
        step_point_bc=reshape(step_point_bc,(1,3))
    else:
        step_point_bc=empty((step2,3))
    
    if all(isnan(step_point_bc))==False:
        for n in range(step2):
            if n==0:
                step_point_bc[n,0]=((1/(2*step2))*Points[tri_single_c,0])+(1-(1/(2*step2)))*Points[tri_single_b,0]
                step_point_bc[n,1]=((1/(2*step2))*Points[tri_single_c,1])+(1-(1/(2*step2)))*Points[tri_single_b,1]
             
                step_point_bc_lastx=((1/(2*step2))*Points[tri_single_b,0])+(1-(1/(2*step2)))*Points[tri_single_c,0]
                step_point_bc_lasty=((1/(2*step2))*Points[tri_single_b,1])+(1-(1/(2*step2)))*Points[tri_single_c,1]
            elif n==(step2-1):
                step_point_bc[n,0]=step_point_bc_lastx
                step_point_bc[n,1]=step_point_bc_lasty    
            else:
                step_point_bc[n,0]=((n/(step2-1))*step_point_bc_lastx)+(1-(n/(step2-1)))*step_point_bc[0,0]
                step_point_bc[n,1]=((n/(step2-1))*step_point_bc_lasty)+(1-(n/(step2-1)))*step_point_bc[0,1]
        
        
        
            if Points[tri_single_b,2]< Points[tri_single_c,2]:
                step_point_bc[n,2]=Points[tri_single_b,2]+(n+1)*step_size
            elif Points[tri_single_c,2]< Points[tri_single_b,2]:
                step_point_bc[n,2]=Points[tri_single_c,2]+(n+1)*step_size
            else:
                step_point_bc[n,2]=Points[tri_single_b,2]
        if Points[tri_single_c,2]< Points[tri_single_b,2]:
            step_point_bc[:,2]=flip(step_point_bc[:,2],0)
        else:
            step_point_bc[:,2]=step_point_bc[:,2]
    else:
        pass
    

#%%  Countour Lines for Triangle Side CA
    step3=int(abs(Points[tri_single_c,2]-Points[tri_single_a,2])/step_size)
    
    if step3==0:
        step_point_ca=array([nan,nan,nan])
        step_point_ca=reshape(step_point_ca,(1,3))
    else:
        step_point_ca=empty((step3,3))
    
    if all(isnan(step_point_ca))==False:
        for n in range(step3):
            if n==0:
                step_point_ca[n,0]=((1/(2*step3))*Points[tri_single_a,0])+(1-(1/(2*step3)))*Points[tri_single_c,0]
                step_point_ca[n,1]=((1/(2*step3))*Points[tri_single_a,1])+(1-(1/(2*step3)))*Points[tri_single_c,1]
             
                step_point_ca_lastx=((1/(2*step3))*Points[tri_single_c,0])+(1-(1/(2*step3)))*Points[tri_single_a,0]
                step_point_ca_lasty=((1/(2*step3))*Points[tri_single_c,1])+(1-(1/(2*step3)))*Points[tri_single_a,1]
            elif n==(step3-1):
                step_point_ca[n,0]=step_point_ca_lastx
                step_point_ca[n,1]=step_point_ca_lasty    
            else:
                step_point_ca[n,0]=((n/(step3-1))*step_point_ca_lastx)+(1-(n/(step3-1)))*step_point_ca[0,0]
                step_point_ca[n,1]=((n/(step3-1))*step_point_ca_lasty)+(1-(n/(step3-1)))*step_point_ca[0,1]
        
            if Points[tri_single_c,2]< Points[tri_single_a,2]:
                step_point_ca[n,2]=Points[tri_single_c,2]+(n+1)*step_size
            elif Points[tri_single_a,2]< Points[tri_single_c,2]:
                step_point_ca[n,2]=Points[tri_single_a,2]+(n+1)*step_size
            else:
                step_point_ca[n,2]=Points[tri_single_c,2]

        if Points[tri_single_a,2]< Points[tri_single_c,2]:
            step_point_ca[:,2]=flip(step_point_ca[:,2],0)
        else:
            step_point_ca[:,2]=step_point_ca[:,2]
    else:
        pass
    
#%%
    if size(step_point_ab,0)>size(step_point_bc,0) and size(step_point_ab,0)>size(step_point_ca,0):
        step_point=step_point_ab
        
        if abs(step_point_bc[step2-1,2]-step_point_ca[step3-1,2])==0 or abs(step_point_bc[step2-1,2]-step_point_ca[step3-1,2])==step_size:
            step_point_opposite=vstack((step_point_bc,step_point_ca))
        elif abs(step_point_bc[step2-1,2]-step_point_ca[0,2])==0 or abs(step_point_bc[step2-1,2]-step_point_ca[0,2])==step_size:
            step_point_opposite=vstack((step_point_bc,step_point_ca))        
        elif abs(step_point_bc[0,2]-step_point_ca[step3-1,2])==0 or abs(step_point_ca[0,2]-step_point_bc[step2-1,2])==step_size:
            step_point_opposite=vstack((step_point_ca,step_point_bc))
        else:
            step_point_opposite=vstack((step_point_ca,flip(step_point_bc,0)))   
    
    elif size(step_point_bc,0)>size(step_point_ab,0) and size(step_point_bc,0)>size(step_point_ca,0):
        step_point=step_point_bc
        
        if abs(step_point_ca[step3-1,2]-step_point_ab[step1-1,2])==0 or abs(step_point_ca[step3-1,2]-step_point_ab[step1-1,2])==step_size:
            step_point_opposite=vstack((step_point_ca,step_point_ab))
        elif abs(step_point_ca[step3-1,2]-step_point_ab[0,2])==0 or abs(step_point_ca[step3-1,2]-step_point_ab[0,2])==step_size:
            step_point_opposite=vstack((step_point_ca,step_point_ab))        
        elif abs(step_point_ca[0,2]-step_point_ab[step1-1,2])==0 or abs(step_point_ca[0,2]-step_point_ab[step1-1,2])==step_size:
            step_point_opposite=vstack((step_point_ab,step_point_ca))
        else:
            step_point_opposite=vstack((step_point_ab,flip(step_point_ca,0)))           
    
    elif size(step_point_ca,0)>size(step_point_ab,0) and size(step_point_ca,0)>size(step_point_bc,0):
        step_point=step_point_ca
        
        if abs(step_point_ab[step1-1,2]-step_point_bc[step2-1,2])==0 or abs(step_point_ab[step1-1,2]-step_point_bc[step2-1,2])==step_size:
            step_point_opposite=vstack((step_point_ab,step_point_bc))
        elif abs(step_point_ab[step1-1,2]-step_point_bc[0,2])==0 or abs(step_point_ab[step1-1,2]-step_point_bc[0,2])==step_size:
            step_point_opposite=vstack((step_point_ab,step_point_bc))        
        elif abs(step_point_ab[0,2]-step_point_bc[step2-1,2])==0 or abs(step_point_bc[0,2]-step_point_ab[step1-1,2])==step_size:
            step_point_opposite=vstack((step_point_bc,step_point_ab))
        else:
            step_point_opposite=vstack((step_point_bc,flip(step_point_ab,0)))   
    
    elif size(step_point_ab,0)==size(step_point_bc,0) and any(isnan(step_point_ca))==True and any(isnan(step_point_ab))==False:
        step_point=step_point_ab
        
        if step_point_ab[0,2]==step_point[0,2]:
            step_point_opposite=step_point_bc
        else:
            step_point_opposite=flip(step_point_bc,0)
    
    elif size(step_point_bc,0)==size(step_point_ca,0) and any(isnan(step_point_ab))==True and any(isnan(step_point_bc))==False:
        step_point=step_point_bc
        
        if step_point_ca[0,2]==step_point[0,2]:
            step_point_opposite=step_point_ca
        else:
            step_point_opposite=flip(step_point_ca,0)    
    
    elif size(step_point_ca,0)==size(step_point_ab,0) and any(isnan(step_point_bc))==True and any(isnan(step_point_ca))==False:
        step_point=step_point_ca
        
        if step_point_ab[0,2]==step_point[0,2]:
            step_point_opposite=step_point_ab
        else:
            step_point_opposite=flip(step_point_ab,0)
    
    else:
        pass
#%%
    if step_point[0,2]==step_point_opposite[0,2]:
        step_point=step_point
        step_point_opposite=step_point_opposite
    else:
        
        step_point=step_point
        step_point_opposite=flip(step_point_opposite,0)
#%%
       
    for n in range(size(step_point,0)):
       ax.plot((step_point[n,0],step_point_opposite[n,0]),(step_point[n,1],step_point_opposite[n,1]),'-k',linewidth=0.1)
    
    Line_values=empty((size(vstack((step_point,step_point_opposite,step_point)),0),2),dtype=int)
 

#    for j in range(3*size(step_point,0)):
#        if any(isin(arange(0,3*size(step_point,0),3),j))==True:
#            for n in range(size(step_point,0)):
#                Line_values[j,:]=str(step_point[n,0:2])
#        elif any(isin(arange(1,3*size(step_point_opposite,0),3),j))==True:
#            for n in range(size(step_point_opposite,0)):
#                Line_values[j,:]=str(step_point_opposite[n,0:2])
#        else:
#            Line_values[j,:]=nan
    
    #Line_values=savetxt("Line_values.scr",Line_values,delimiter=",",header="_multiple _line",comments="")

    
    
    
    tri_singles_a=tri_single_a+((empty((size(triangles,0),1),dtype=int))*0)
    tri_singles_b=tri_single_b+(empty((size(triangles,0),1),dtype=int)*0)
    tri_singles_c=tri_single_c+(empty((size(triangles,0),1),dtype=int)*0)
    
    tri_single=hstack((tri_singles_a,tri_singles_b,tri_singles_c))
    
    #ax.triplot(Points[:,0],Points[:,1],tri_single,'--b',linewidth=0.05)
    
    #plt.pause(0.15)

plt.show()

plt.savefig("Site1.pdf",transparent=True)


