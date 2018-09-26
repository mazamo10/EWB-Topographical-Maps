#%% Site 1 Topographical Map--Matthew Zamora

#%%

from numpy import *
import pandas as pd
from scipy.optimize import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.spatial import Delaunay

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
B=math.acos((b**2-a**2-c**2)/(-2*a*c))            # Angle in radians opposite b

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

survey=pd.read_excel('Topo1.xlsx')
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
# A new theta is calculated by adjusting the axis so that the x axis lies
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



#%% 2D Point Plot

fig1=plt.figure()
ax = fig1.add_subplot(111)


for i in range(size(X)):
    if Key_Feature[i,0]==1:
        ax.scatter(X[i],Y[i],c='b')
    else:
        ax.scatter(X[i],Y[i],c='k')
for x, y, z in zip(X, Y, Z):
    # Annotate the points 5 _points_ above and to the left of the vertex
    ax.annotate('{}'.format(z), xy=(x,y), xytext=(-5, 5), ha='right',
                textcoords='offset points')


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


step_size=1/2;

step_length=empty((size(triangles,0),1))

   
for i in range(size(triangles,0)):
    tri_a=int(triangles[i,0])
    tri_b=int(triangles[i,1])
    tri_c=int(triangles[i,2])    
    
    step1=int(abs(Points[tri_a,2]-Points[tri_b,2])/step_size)    
    step2=int(abs(Points[tri_b,2]-Points[tri_c,2])/step_size)
    step3=int(abs(Points[tri_c,2]-Points[tri_a,2])/step_size)
    
    step_length[i,0]=max((step1,step2,step3))

max_length=int(max(step_length))

contours=abs(empty((size(triangles,0),max_length,4))*0)

for i in range(size(triangles,0)):
    tri_a=int(triangles[i,0])
    tri_b=int(triangles[i,1])
    tri_c=int(triangles[i,2])
    

#%%  Countour Lines for Triangle Side AB  
    step1=int(abs(Points[tri_a,2]-Points[tri_b,2])/step_size)
    
    if step1==0:
        step_ab=array([nan,nan,nan])
        step_ab=reshape(step_ab,(1,3))
    else:
        step_ab=empty((step1,3))
    if all(isnan(step_ab))==False:
        for n in range(step1):
            if n==0:
                step_ab[n,0]=((1/(2*step1))*Points[tri_b,0])+(1-(1/(2*step1)))*Points[tri_a,0] # Creates value closest to Points[tri_a,0]
                step_ab[n,1]=((1/(2*step1))*Points[tri_b,1])+(1-(1/(2*step1)))*Points[tri_a,1] # Creates value closest to Points[tri_a,1]
             
                step_point_ab_lastx=((1/(2*step1))*Points[tri_a,0])+(1-(1/(2*step1)))*Points[tri_b,0]# Creates value closest to Points[tri_b,0]
                step_point_ab_lasty=((1/(2*step1))*Points[tri_a,1])+(1-(1/(2*step1)))*Points[tri_b,1]# Creates value closest to Points[tri_b,1]
             
            elif n==(step1-1):
                step_ab[n,0]=step_point_ab_lastx # Assigns x value of last point
                step_ab[n,1]=step_point_ab_lasty # Assigns y value of last point    
            else:
                step_ab[n,0]=((n/(step1-1))*step_point_ab_lastx)+(1-(n/(step1-1)))*step_ab[0,0] # Determines x values of step points between the fisrst and last step point value iterating from points closest to the first value
                step_ab[n,1]=((n/(step1-1))*step_point_ab_lasty)+(1-(n/(step1-1)))*step_ab[0,1] # Determines y values of step points between the fisrst and last step point value iterating from points closest to the first value
    
            if Points[tri_a,2]< Points[tri_b,2]: # If point a is lower than point b
                step_ab[n,2]=Points[tri_a,2]+(n+1)*step_size # Creates z values that increase between a and b
            elif Points[tri_b,2]< Points[tri_a,2]: # If point b is lower than point a
                step_ab[n,2]=Points[tri_b,2]+(n+1)*step_size # Creates z values that decrease between a and b
            else: # If z values of points a and b are equal
                step_ab[n,2]=Points[tri_b,2] # z value is the same as both points
        if Points[tri_b,2]< Points[tri_a,2]:
            step_ab[:,2]=flip(step_ab[:,2],0)
        else:
            step_ab[:,2]=step_ab[:,2]
    else:
         pass
     
#%%  Countour Lines for Triangle Side BC  
    step2=int(abs(Points[tri_b,2]-Points[tri_c,2])/step_size)
    
    if step2==0:
        step_bc=array([nan,nan,nan])
        step_bc=reshape(step_bc,(1,3))
    else:
        step_bc=empty((step2,3))
    
    if all(isnan(step_bc))==False:
        for n in range(step2):
            if n==0:
                step_bc[n,0]=((1/(2*step2))*Points[tri_c,0])+(1-(1/(2*step2)))*Points[tri_b,0]
                step_bc[n,1]=((1/(2*step2))*Points[tri_c,1])+(1-(1/(2*step2)))*Points[tri_b,1]
             
                step_point_bc_lastx=((1/(2*step2))*Points[tri_b,0])+(1-(1/(2*step2)))*Points[tri_c,0]
                step_point_bc_lasty=((1/(2*step2))*Points[tri_b,1])+(1-(1/(2*step2)))*Points[tri_c,1]
            elif n==(step2-1):
                step_bc[n,0]=step_point_bc_lastx
                step_bc[n,1]=step_point_bc_lasty    
            else:
                step_bc[n,0]=((n/(step2-1))*step_point_bc_lastx)+(1-(n/(step2-1)))*step_bc[0,0]
                step_bc[n,1]=((n/(step2-1))*step_point_bc_lasty)+(1-(n/(step2-1)))*step_bc[0,1]
        
        
        
            if Points[tri_b,2]< Points[tri_c,2]:
                step_bc[n,2]=Points[tri_b,2]+(n+1)*step_size
            elif Points[tri_c,2]< Points[tri_b,2]:
                step_bc[n,2]=Points[tri_c,2]+(n+1)*step_size
            else:
                step_bc[n,2]=Points[tri_b,2]
        if Points[tri_c,2]< Points[tri_b,2]:
            step_bc[:,2]=flip(step_bc[:,2],0)
        else:
            step_bc[:,2]=step_bc[:,2]
    else:
        pass
    

#%%  Countour Lines for Triangle Side CA
    step3=int(abs(Points[tri_c,2]-Points[tri_a,2])/step_size)
    
    if step3==0:
        step_ca=array([nan,nan,nan])
        step_ca=reshape(step_ca,(1,3))
    else:
        step_ca=empty((step3,3))
    
    if all(isnan(step_ca))==False:
        for n in range(step3):
            if n==0:
                step_ca[n,0]=((1/(2*step3))*Points[tri_a,0])+(1-(1/(2*step3)))*Points[tri_c,0]
                step_ca[n,1]=((1/(2*step3))*Points[tri_a,1])+(1-(1/(2*step3)))*Points[tri_c,1]
             
                step_point_ca_lastx=((1/(2*step3))*Points[tri_c,0])+(1-(1/(2*step3)))*Points[tri_a,0]
                step_point_ca_lasty=((1/(2*step3))*Points[tri_c,1])+(1-(1/(2*step3)))*Points[tri_a,1]
            elif n==(step3-1):
                step_ca[n,0]=step_point_ca_lastx
                step_ca[n,1]=step_point_ca_lasty    
            else:
                step_ca[n,0]=((n/(step3-1))*step_point_ca_lastx)+(1-(n/(step3-1)))*step_ca[0,0]
                step_ca[n,1]=((n/(step3-1))*step_point_ca_lasty)+(1-(n/(step3-1)))*step_ca[0,1]
        
            if Points[tri_c,2]< Points[tri_a,2]:
                step_ca[n,2]=Points[tri_c,2]+(n+1)*step_size
            elif Points[tri_a,2]< Points[tri_c,2]:
                step_ca[n,2]=Points[tri_a,2]+(n+1)*step_size
            else:
                step_ca[n,2]=Points[tri_c,2]

        if Points[tri_a,2]< Points[tri_c,2]:
            step_ca[:,2]=flip(step_ca[:,2],0)
        else:
            step_ca[:,2]=step_ca[:,2]
    else:
        pass
    
#%%
    if size(step_ab,0)>size(step_bc,0) and size(step_ab,0)>size(step_ca,0):
        step_point=step_ab
        
        if abs(step_bc[step2-1,2]-step_ca[step3-1,2])==0 or abs(step_bc[step2-1,2]-step_ca[step3-1,2])==step_size:
            step_point_opposite=vstack((step_bc,step_ca))
        elif abs(step_bc[step2-1,2]-step_ca[0,2])==0 or abs(step_bc[step2-1,2]-step_ca[0,2])==step_size:
            step_point_opposite=vstack((step_bc,step_ca))        
        elif abs(step_bc[0,2]-step_ca[step3-1,2])==0 or abs(step_ca[0,2]-step_bc[step2-1,2])==step_size:
            step_point_opposite=vstack((step_ca,step_bc))
        else:
            step_point_opposite=vstack((step_ca,flip(step_bc,0)))   
    
    elif size(step_bc,0)>size(step_ab,0) and size(step_bc,0)>size(step_ca,0):
        step_point=step_bc
        
        if abs(step_ca[step3-1,2]-step_ab[step1-1,2])==0 or abs(step_ca[step3-1,2]-step_ab[step1-1,2])==step_size:
            step_point_opposite=vstack((step_ca,step_ab))
        elif abs(step_ca[step3-1,2]-step_ab[0,2])==0 or abs(step_ca[step3-1,2]-step_ab[0,2])==step_size:
            step_point_opposite=vstack((step_ca,step_ab))        
        elif abs(step_ca[0,2]-step_ab[step1-1,2])==0 or abs(step_ca[0,2]-step_ab[step1-1,2])==step_size:
            step_point_opposite=vstack((step_ab,step_ca))
        else:
            step_point_opposite=vstack((step_ab,flip(step_ca,0)))           
    
    elif size(step_ca,0)>size(step_ab,0) and size(step_ca,0)>size(step_bc,0):
        step_point=step_ca
        
        if abs(step_ab[step1-1,2]-step_bc[step2-1,2])==0 or abs(step_ab[step1-1,2]-step_bc[step2-1,2])==step_size:
            step_point_opposite=vstack((step_ab,step_bc))
        elif abs(step_ab[step1-1,2]-step_bc[0,2])==0 or abs(step_ab[step1-1,2]-step_bc[0,2])==step_size:
            step_point_opposite=vstack((step_ab,step_bc))        
        elif abs(step_ab[0,2]-step_bc[step2-1,2])==0 or abs(step_bc[0,2]-step_ab[step1-1,2])==step_size:
            step_point_opposite=vstack((step_bc,step_ab))
        else:
            step_point_opposite=vstack((step_bc,flip(step_ab,0)))   
    
    elif size(step_ab,0)==size(step_bc,0) and any(isnan(step_ca))==True and any(isnan(step_ab))==False:
        step_point=step_ab
        
        if step_ab[0,2]==step_point[0,2]:
            step_point_opposite=step_bc
        else:
            step_point_opposite=flip(step_bc,0)
    
    elif size(step_bc,0)==size(step_ca,0) and any(isnan(step_ab))==True and any(isnan(step_bc))==False:
        step_point=step_bc
        
        if step_ca[0,2]==step_point[0,2]:
            step_point_opposite=step_ca
        else:
            step_point_opposite=flip(step_ca,0)    
    
    elif size(step_ca,0)==size(step_ab,0) and any(isnan(step_bc))==True and any(isnan(step_ca))==False:
        step_point=step_ca
        
        if step_ab[0,2]==step_point[0,2]:
            step_point_opposite=step_ab
        else:
            step_point_opposite=flip(step_ab,0)
    
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
       
       

       contours[i,n,0]=step_point[n,0]
       contours[i,n,1]=step_point[n,1]
       contours[i,n,2]=step_point_opposite[n,0]
       contours[i,n,3]=step_point_opposite[n,1]
       

ax.plot(Points[:,0], Points[:,1],'ok',markersize=1)
plt.show()

fig3.savefig("Contour_map.pdf", bbox_inches='tight')
