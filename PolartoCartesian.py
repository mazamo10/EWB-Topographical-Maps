def PolartoCartesian(filename,sheetnumber):
    """
    Created on 10/29/18
    
    @author: Ryan Stephenson
    @author: Matthew Zamora"""
    

    import numpy as np
    from pandas import read_excel


    instrument_height=5.1;
     
    #Unknown Variables
    hypotenuse=[]
    height = []
    angle = []
    
    distance=[]


    x=[]
    y=[]
    z=[]
        
    def pol2cart(distance, angle):
        """Converts polar data to xy coordinates for plotting xyz"""
        r = np.array(distance)
        theta = np.deg2rad(360 - np.array(angle))
         
        x = r * np.cos(theta)
        y = r * np.sin(theta)
         
        return(x, y)

    for i in range(sheetnumber):
        X=[]
        Y=[]
        Z=[]

        survey=read_excel(filename,sheet_name=i)
        survey.values

        hypotenuse=survey['hypotenuse']
        hypotenuse=hypotenuse.values
        
        height=survey['height']
        height=height.values
     
        angle=survey['angle']
        angle=angle.values
            
        distance=np.sqrt(np.power(hypotenuse,2)-np.power(height,2))
        
        dx, dy = pol2cart(distance,angle) 
        dz = -(height)

        if np.size(x)==0:        
            X=dx
            Y=dy
            Z=dz
        else:
            X=x[-1]+dx
            Y=y[-1]+dy
            Z=z[-1]+instrument_height+dz


        x.extend(X)
        y.extend(Y)
        z.extend(Z)   
        
    x=np.reshape(x,(np.size(x),1))
    y=np.reshape(y,(np.size(y),1))
    z=np.reshape(z,(np.size(z),1))
     
    Points=np.hstack((x,y,z))
    return Points
