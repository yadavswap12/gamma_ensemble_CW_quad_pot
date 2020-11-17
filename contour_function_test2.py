# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 17:15:47 2019

@author: yadav
"""
"""
Created on Sat Jun 22 10:52:20 2019

@author: yadav
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 17:09:12 2018

@author: yadav
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 17:09:12 2018

@author: yadav
"""

# code to generate contour nu, nu1 and nu2 for J(s)=(c1s+c0)((s+1)/s)**(1/theta)
import math
import cmath
import numpy

#rho = 2.0    # V(x)=rho*x
#c = theta/rho    # eqn.(4.21) in C.R.
c1=2.0
c0=1.0

n = 10000000    # number of points in each quadrant
    
#    y1 = numpy.linspace(10.0/n,10.0,n)    #linspace creates array of x values, y_max=1.05 for c1=1,c0=1/2  
f_out=file("input_output_files/rho_2.0/contour/nu_contour_function_soft_edge_18000points.txt","w")
#    f_out=open("input_output_files//contour//nu_contour_function_soft_edge_18000points.txt","w")    # open() can also be used to create file    
y1 = numpy.linspace(2.0/n,2.0,n)    #linspace creates array of x values, y_max=1.05 for c1=1,c0=1/2  

for i in range(len(y1)):    # function len() on array gives no. of rows of array
    yy = y1[i]
    if ((1.0/4)+(yy/math.tan(c1*yy))-yy**2)<=0.0:
        r = yy+ 2.0/n
        r = yy
        print r
        break

for i in range(1,18000,1):    #only integer step argument is available for range() hence 0 to 36000. to get 180 degrees will later divide by 200
    y = r*math.sin(math.radians(i/100.0))
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i<=9000:
        x = cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative
        f_out.write(str(x)+" "+str(y)+'\n')
#            f_out1.write(str(x)+" "+str(y)+'\n')
    
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i>9000:
        x = -cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative    
        f_out.write(str(x)+" "+str(y)+'\n')
#            f_out1.write(str(x)+" "+str(y)+'\n')
#f_out.close()    # () at the end is necessary to close the file 
#f_out1.close()    # () at the end is necessary to close the file
#y2 = numpy.linspace(0.7,0,n)
 
#y2 = numpy.linspace(0.062,0,n)     
for i in range(17999,0,-1):    #only integer step argument is available for range() hence 0 to 36000. to get 180 degrees will later divide by 200
    y = -r*math.sin(math.radians(i/100.0))
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i<=9000:
        x = cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative
        f_out.write(str(x)+" "+str(y)+'\n')
#            f_out2.write(str(x)+" "+str(y)+'\n')
    
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i>9000:
        x = -cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative    
        f_out.write(str(x)+" "+str(y)+'\n')
#            f_out2.write(str(x)+" "+str(y)+'\n')

        
f_out.close()    # () at the end is necessary to close the file
