#%% Site 2 Topographical Map--Matthew Zamora

#%%

from numpy import *
import matplotlib.pyplot as plt
from contours import *

#%% Refrence values
a=40.5                              # Distance in feet between A and instrument
b=68                                # Distance in feet between B and instrument
theta_ref_A=139                     # Reference angle in degrees measured by instrument to reference point A
theta_ref_B=169                     # Reference angle in degrees measured by instrument to reference point B
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

survey=pd.read_excel('Topo2.xlsx')
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
#%%

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
Contours(Points,X,Y,Z,1/2)