#!/usr/bin/env python3
# -*- coding: utf-8 -*-
def Contours(x,y,z,step_size,file_saved):

    """
    Created on Spring 2019
    
    @author: Mattew Zamora
    
    Admittedly the code for this function is gross and overcomplicated; however, it still manages to complete its assigned task.
    If anyone in later years is looking back on this, I'm sorry for the mess, please feel free to simplify it as much as possible.
    
    The x,y,z values of the survey points are passed to this function along with the desired step size (at its current state the function
    tends to break down for step size values greater than 1). The function then performs a common triangulation algorithm on the points which
    creates an array of indices denoting which set of three points combined to form each triangle. These triangles are then iterated through.
    For each iteration, the height difference between each endpoint of the triangle is calculated and additional points are linearly interpolated
    so that two adjacent points along the same line signify a height difference equal to one step size. This creates three separate arrays
    (step_ab, step_bc, step_ca) containing the x,y,z coordinates of the interpolated points. The rest of the code is devoted to organizing these
    arrays into two arrays of equal length (in most cases) where each successive value corresponds to the same height for both arrays. Finally,
    the function returns a plot of the produced map and a pdf file of the map.
    """

# Imports
    from scipy.spatial import Delaunay
    import numpy as np
    import matplotlib.pyplot as plt
# Triangulation   
    points=np.hstack((x,y))
    tri=Delaunay(points)

    triangles=tri.simplices
    points=np.hstack((x,y,z))

    fig=plt.figure()
    plt.gca().set_aspect('equal',adjustable='box')
    ax = fig.add_subplot(111)


    for i in range(np.size(triangles,0)):
        tri_a=int(triangles[i,0]);tri_b=int(triangles[i,1]);tri_c=int(triangles[i,2])
    
# Countour Lines for Triangle Side AB  
        step1=int(abs(points[tri_a,2]-points[tri_b,2])/step_size)
    
        if step1==0:
            step_ab=np.array([np.nan,np.nan,np.nan])
            step_ab=np.reshape(step_ab,(1,3))
        else:
            step_ab=np.empty((step1,3))
        
            for n in range(step1):
                if n==0:
                    step_ab[n,0]=((1/(2*step1))*points[tri_b,0])+(1-(1/(2*step1)))*points[tri_a,0]
                    step_ab[n,1]=((1/(2*step1))*points[tri_b,1])+(1-(1/(2*step1)))*points[tri_a,1]
             
                    step_point_ab_lastx=((1/(2*step1))*points[tri_a,0])+(1-(1/(2*step1)))*points[tri_b,0]
                    step_point_ab_lasty=((1/(2*step1))*points[tri_a,1])+(1-(1/(2*step1)))*points[tri_b,1]
                 
                elif n==(step1-1):
                    step_ab[n,0]=step_point_ab_lastx
                    step_ab[n,1]=step_point_ab_lasty  
                else:
                    step_ab[n,0]=((n/(step1-1))*step_point_ab_lastx)+(1-(n/(step1-1)))*step_ab[0,0]
                    step_ab[n,1]=((n/(step1-1))*step_point_ab_lasty)+(1-(n/(step1-1)))*step_ab[0,1]
        
                if points[tri_a,2]< points[tri_b,2]:
                    step_ab[n,2]=points[tri_a,2]+(n+1)*step_size
                elif points[tri_b,2]< points[tri_a,2]:
                    step_ab[n,2]=points[tri_b,2]+(n+1)*step_size
                else:
                    step_ab[n,2]=points[tri_b,2]
            if points[tri_b,2]< points[tri_a,2]:
                step_ab[:,2]=np.flip(step_ab[:,2],0)
            else:
                step_ab[:,2]=step_ab[:,2]      
         
# Countour Lines for Triangle Side BC  
        step2=int(abs(points[tri_b,2]-points[tri_c,2])/step_size)
        
        if step2==0:
            step_bc=np.array([np.nan,np.nan,np.nan])
            step_bc=np.reshape(step_bc,(1,3))
        else:
            step_bc=np.empty((step2,3))
        
            for n in range(step2):
                if n==0:
                    step_bc[n,0]=((1/(2*step2))*points[tri_c,0])+(1-(1/(2*step2)))*points[tri_b,0]
                    step_bc[n,1]=((1/(2*step2))*points[tri_c,1])+(1-(1/(2*step2)))*points[tri_b,1]
                 
                    step_point_bc_lastx=((1/(2*step2))*points[tri_b,0])+(1-(1/(2*step2)))*points[tri_c,0]
                    step_point_bc_lasty=((1/(2*step2))*points[tri_b,1])+(1-(1/(2*step2)))*points[tri_c,1]
                elif n==(step2-1):
                    step_bc[n,0]=step_point_bc_lastx
                    step_bc[n,1]=step_point_bc_lasty    
                else:
                    step_bc[n,0]=((n/(step2-1))*step_point_bc_lastx)+(1-(n/(step2-1)))*step_bc[0,0]
                    step_bc[n,1]=((n/(step2-1))*step_point_bc_lasty)+(1-(n/(step2-1)))*step_bc[0,1]
            
            
            
                if points[tri_b,2]< points[tri_c,2]:
                    step_bc[n,2]=points[tri_b,2]+(n+1)*step_size
                elif points[tri_c,2]< points[tri_b,2]:
                    step_bc[n,2]=points[tri_c,2]+(n+1)*step_size
                else:
                    step_bc[n,2]=points[tri_b,2]
            if points[tri_c,2]< points[tri_b,2]:
                step_bc[:,2]=np.flip(step_bc[:,2],0)
            else:
                step_bc[:,2]=step_bc[:,2]
        
# Countour Lines for Triangle Side CA
        step3=int(abs(points[tri_c,2]-points[tri_a,2])/step_size)
        
        if step3==0:
            step_ca=np.array([np.nan,np.nan,np.nan])
            step_ca=np.reshape(step_ca,(1,3))
        else:
            step_ca=np.empty((step3,3))
        
            for n in range(step3):
                if n==0:
                    step_ca[n,0]=((1/(2*step3))*points[tri_a,0])+(1-(1/(2*step3)))*points[tri_c,0]
                    step_ca[n,1]=((1/(2*step3))*points[tri_a,1])+(1-(1/(2*step3)))*points[tri_c,1]
                 
                    step_point_ca_lastx=((1/(2*step3))*points[tri_c,0])+(1-(1/(2*step3)))*points[tri_a,0]
                    step_point_ca_lasty=((1/(2*step3))*points[tri_c,1])+(1-(1/(2*step3)))*points[tri_a,1]
                elif n==(step3-1):
                    step_ca[n,0]=step_point_ca_lastx
                    step_ca[n,1]=step_point_ca_lasty    
                else:
                    step_ca[n,0]=((n/(step3-1))*step_point_ca_lastx)+(1-(n/(step3-1)))*step_ca[0,0]
                    step_ca[n,1]=((n/(step3-1))*step_point_ca_lasty)+(1-(n/(step3-1)))*step_ca[0,1]
            
                if points[tri_c,2]< points[tri_a,2]:
                    step_ca[n,2]=points[tri_c,2]+(n+1)*step_size
                elif points[tri_a,2]< points[tri_c,2]:
                    step_ca[n,2]=points[tri_a,2]+(n+1)*step_size
                else:
                    step_ca[n,2]=points[tri_c,2]
    
            if points[tri_a,2]< points[tri_c,2]:
                step_ca[:,2]=np.flip(step_ca[:,2],0)
            else:
                step_ca[:,2]=step_ca[:,2]
        
# Sorting points into appropriate order
        if np.size(step_ab,0)>np.size(step_bc,0) and np.size(step_ab,0)>np.size(step_ca,0):
            step_point=step_ab
            
            if step_ab[0,2]==step_bc[0,2] and step_ab[step1-1,2]==step_ca[step3-1,2]:
                step_point_opposite=np.vstack((step_bc,step_ca))
            elif step_ab[0,2]==step_bc[0,2] and step_ab[step1-1,2]==step_ca[0,2]:
                step_point_opposite=np.vstack((step_bc,np.flip(step_ca,0)))        
            elif step_ab[0,2]==step_bc[step2-1,2] and step_ab[step1-1,2]==step_ca[step3-1,2]:
                step_point_opposite=np.vstack((np.flip(step_bc,0),step_ca))
            elif step_ab[0,2]==step_bc[step2-1,2] and step_ab[step1-1,2]==step_ca[0,2]:
                step_point_opposite=np.vstack((np.flip(step_bc,0),np.flip(step_ca,0)))
            
            elif step_ab[0,2]==step_ca[0,2] and step_ab[step1-1,2]==step_bc[step2-1,2]:
                step_point_opposite=np.vstack((step_ca,step_bc))
            elif step_ab[0,2]==step_ca[0,2] and step_ab[step1-1,2]==step_bc[0,2]:
                step_point_opposite=np.vstack((step_ca,np.flip(step_bc,0)))        
            elif step_ab[0,2]==step_ca[step3-1,2] and step_ab[step1-1,2]==step_bc[step2-1,2]:
                step_point_opposite=np.vstack((np.flip(step_ca,0),step_bc))
            else:
                step_point_opposite=np.vstack((np.flip(step_ca,0),np.flip(step_bc,0)))            

        
        elif np.size(step_bc,0)>np.size(step_ab,0) and np.size(step_bc,0)>np.size(step_ca,0):
            step_point=step_bc
            
            if step_bc[0,2]==step_ab[0,2] and step_bc[step2-1,2]==step_ca[step3-1,2]:
                step_point_opposite=np.vstack((step_ab,step_ca))
            elif step_bc[0,2]==step_ab[0,2] and step_bc[step2-1,2]==step_ca[0,2]:
                step_point_opposite=np.vstack((step_ab,np.flip(step_ca,0)))        
            elif step_bc[0,2]==step_ab[step1-1,2] and step_bc[step2-1,2]==step_ca[step3-1,2]:
                step_point_opposite=np.vstack((np.flip(step_ab,0),step_ca))
            elif step_bc[0,2]==step_ab[step1-1,2] and step_bc[step2-1,2]==step_ca[0,2]:
                step_point_opposite=np.vstack((np.flip(step_ab,0),np.flip(step_ca,0)))
            
            elif step_bc[0,2]==step_ca[0,2] and step_bc[step2-1,2]==step_ab[step1-1,2]:
                step_point_opposite=np.vstack((step_ca,step_ab))
            elif step_bc[0,2]==step_ca[0,2] and step_bc[step2-1,2]==step_ab[0,2]:
                step_point_opposite=np.vstack((step_ca,np.flip(step_ab,0)))        
            elif step_bc[0,2]==step_ca[step3-1,2] and step_bc[step2-1,2]==step_ab[step1-1,2]:
                step_point_opposite=np.vstack((np.flip(step_ca,0),step_ab))
            else:
                step_point_opposite=np.vstack((np.flip(step_ca,0),np.flip(step_ab,0)))            

         
        elif np.size(step_ca,0)>np.size(step_ab,0) and np.size(step_ca,0)>np.size(step_bc,0):
            step_point=step_ca
            
            if step_ca[0,2]==step_ab[0,2] and step_ca[step3-1,2]==step_bc[step2-1,2]:
                step_point_opposite=np.vstack((step_ab,step_bc))
            elif step_ca[0,2]==step_ab[0,2] and step_ca[step3-1,2]==step_bc[0,2]:
                step_point_opposite=np.vstack((step_ab,np.flip(step_bc,0)))        
            elif step_ca[0,2]==step_ab[step1-1,2] and step_ca[step3-1,2]==step_bc[step2-1,2]:
                step_point_opposite=np.vstack((np.flip(step_ab,0),step_bc))
            elif step_ca[0,2]==step_ab[step1-1,2] and step_ca[step3-1,2]==step_bc[0,2]:
                step_point_opposite=np.vstack((np.flip(step_ab,0),np.flip(step_bc,0)))
            
            elif step_ca[0,2]==step_bc[0,2] and step_ca[step3-1,2]==step_ab[step1-1,2]:
                step_point_opposite=np.vstack((step_bc,step_ab))
            elif step_ca[0,2]==step_bc[0,2] and step_ca[step3-1,2]==step_ab[0,2]:
                step_point_opposite=np.vstack((step_bc,np.flip(step_ab,0)))        
            elif step_ca[0,2]==step_bc[step2-1,2] and step_ca[step3-1,2]==step_ab[step1-1,2]:
                step_point_opposite=np.vstack((np.flip(step_bc,0),step_ab))
            else:
                step_point_opposite=np.vstack((np.flip(step_bc,0),np.flip(step_ab,0))) 
        
        elif np.size(step_ab,0)==np.size(step_bc,0) and np.any(np.isnan(step_ca))==True and np.any(np.isnan(step_ab))==False:
            step_point=step_ab
            
            if step_ab[0,2]==step_point[0,2]:
                step_point_opposite=step_bc
            else:
                step_point_opposite=np.flip(step_bc,0)
        
        elif np.size(step_bc,0)==np.size(step_ca,0) and np.any(np.isnan(step_ab))==True and np.any(np.isnan(step_bc))==False:
            step_point=step_bc
            
            if step_ca[0,2]==step_point[0,2]:
                step_point_opposite=step_ca
            else:
                step_point_opposite=np.flip(step_ca,0)    
        
        elif np.size(step_ca,0)==np.size(step_ab,0) and np.any(np.isnan(step_bc))==True and np.any(np.isnan(step_ca))==False:
            step_point=step_ca
            
            if step_ab[0,2]==step_point[0,2]:
                step_point_opposite=step_ab
            else:
                step_point_opposite=np.flip(step_ab,0)
        else:
            step_point=np.nan;step_point_opposite=np.nan


        if np.any(np.isnan(step_point))==True:
            pass
        elif step_point[0,2]==step_point_opposite[0,2]:
            pass
        else:
           step_point_opposite=np.flip(step_point_opposite,0)
        
        if np.any(np.isnan(step_point))==False:
            for n in range(np.size(step_point,0)):
                ax.plot((step_point[n,0],step_point_opposite[n,0]),(step_point[n,1],step_point_opposite[n,1]),'-k',linewidth=0.1)
        else:
            pass
    ax.plot(points[:,0], points[:,1],'ok',markersize=1)
    plt.show()
    
    fig.savefig(file_saved, bbox_inches='tight')
    return
