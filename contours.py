def Contours(Points,X,Y,Z,step_size):
    
    
    from scipy.spatial import Delaunay
    from numpy import hstack, vstack, size, empty, isnan, reshape, all, array, nan, flip, any
    import pandas as pd
    import matplotlib.pyplot as plt
    
    
    tri=Delaunay(Points)

    triangles=tri.simplices
    Points=hstack((X,Y,Z))

    fig3=plt.figure()
    plt.gca().set_aspect('equal',adjustable='box')
    ax = fig3.add_subplot(111)


    step_length=empty((size(triangles,0),1))

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

    ax.plot(Points[:,0], Points[:,1],'ok',markersize=1)
    plt.show()
    
    fig3.savefig("Contour_map.pdf", bbox_inches='tight')
    return
