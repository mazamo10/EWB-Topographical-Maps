
#Ryan Stephenson 10/29/18
import numpy as np

#Unknown Variables
distance = []
angle = []
height = []
elevation = 0.0

#import data (W.I.P.)
			#  height, angle, distance from Monitoring Site A excel data
data_str_1 = '''5001.6	180	50
				4999.6	331	30.1
				4999.8	312	29.4
				4999.4	296	31.6
				4999.7	275	29.5
				4999.5	261	29.1
				4999.9	249	29.6
				5000.4	233	28.1
				5000.4	208	29
				5000.9	188	31.1
				4999.8	357	62.2
				4999.3	347	59.9
				4999.1	330	59.8
				4999.5	306	58.4
				4999.9	285	58
				5000.3	260	58.5
				5000.8	236	57.9
				5001.2	212	61.4
				5001.2	196	59.6
				4999	352	88.7
				4999.6	318	90.5
				5000.4	292	89.5
				5001.2	253	89
				5001.6	220	90.6
				5001.9	194	88.4'''

data_set_1 = np.fromstring(data_str_1,dtype=float,sep=' ').reshape(25,3) #creates data array
height, angle, distance = np.hsplit(data_set_1,3 ) #format data array


#Converts polar data to xy coordinates for plotting xyz
def pol2cart(distance, angle):
	r = np.array(distance)
	theta = np.deg2rad(360 - np.array(angle))

	x = r * np.cos(theta)
	y = r * np.sin(theta) 

	return(x, y)

x, y = pol2cart(distance,angle)
height = (height - 5000) + elevation

#print height, x, y 