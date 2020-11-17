# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:14:32 2018

@author: yadav
"""

# code to generate contour nu, nu1 and nu2 for J(s)=(c1s+c0)((s+1)/s)**(1/theta)
import math
import cmath
import numpy

#theta = 2.0
#alpha=8.0
rho = 2.0    # V(x)=gamma*(rho*x+alpha*x**(1.0/alpha)), alpha>1
gamma = 0.9
iteration = 8.0

#c = theta/rho    # eqn.(4.21) in C.R.
c1 = 0.94944487275685974    # renormalized c1 for gamma=0.5 computed from python file 'renormalized_joukowsky_parameters_c1_c0.py'
c0 = 0.45183538883128627    # renormalized c0 for gamma=0.5 computed from python file 'renormalized_joukowsky_parameters_c1_c0.py' 

n = 100000    # number of points in each quadrant

f_out=file("input_output_files/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu_contour_gamma="+str(gamma)+"_soft_edge_18000points_iter"+str(iteration)+".txt","w")
f_out1=file("input_output_files/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu1_contour_gamma="+str(gamma)+"_soft_edge_18000points_iter"+str(iteration)+".txt","w")
f_out2=file("input_output_files/rho_"+str(rho)+"/gamma_"+str(gamma)+"/contour/nu2_contour_gamma="+str(gamma)+"_soft_edge_18000points_iter"+str(iteration)+".txt","w")

y1 = numpy.linspace(20.0/n,20.0,n)    #linspace creates array of x values, y_max=1.05 for c1=1,c0=1/2  

for i in range(len(y1)):    # function len() on array gives no. of rows of array
    yy = y1[i]
    if ((1.0/4)+(yy/math.tan(c1*yy))-yy**2)<0.0:
        yy1 = y1[i-1]
        yy2 = y1[i] 
#            r = yy+ 2.0/n
#            yy_avg = (y1[i]+y1[i-1])/2.0
#            if ((1.0/4)+(yy_avg/math.tan(c1*yy_avg))-yy_avg**2)<=0.0:
        while True:
            yy3 = (yy1+yy2)/2.0
            if ((1.0/4)+(yy3/math.tan(c1*yy3))-yy3**2)<0.0:
                yy2 = yy3
                    
            if ((1.0/4)+(yy3/math.tan(c1*yy3))-yy3**2)>=0.0:
                yy1 = yy3
                
            r=yy2
            error = yy1-yy2
            if abs(error)<=(1e-15):
                break
        break


for i in range(1,18000,1):    #only integer step argument is available for range() hence 0 to 36000. to get 180 degrees will later divide by 200
    y = r*math.sin(math.radians(i/100.0))
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i<=9000:
        x = cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative
        f_out.write(str(x)+" "+str(y)+'\n')
        f_out1.write(str(x)+" "+str(y)+'\n')
    
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i>9000:
        x = -cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative    
        f_out.write(str(x)+" "+str(y)+'\n')
        f_out1.write(str(x)+" "+str(y)+'\n')
#f_out.close()    # () at the end is necessary to close the file 
#f_out1.close()    # () at the end is necessary to close the file
#y2 = numpy.linspace(0.7,0,n)
 
#y2 = numpy.linspace(0.062,0,n)     
for i in range(17999,0,-1):    #only integer step argument is available for range() hence 0 to 36000. to get 180 degrees will later divide by 200
    y = -r*math.sin(math.radians(i/100.0))
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i<=9000:
        x = cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative
        f_out.write(str(x)+" "+str(y)+'\n')
        f_out2.write(str(x)+" "+str(y)+'\n')
    
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i>9000:
        x = -cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative    
        f_out.write(str(x)+" "+str(y)+'\n')
        f_out2.write(str(x)+" "+str(y)+'\n')
#f_out.close()    # () at the end is necessary to close the file 
#f_out1.close()    # () at the end is necessary to close the file
#y3 = numpy.linspace(0,-0.7,n)
#y3 = numpy.linspace(0,-0.062,n)      
#for i in range(len(y3)):    
#    y = y3[i]
#    x = (-1)*cmath.sqrt((1.0/4)+(y3[i]/math.tan(c1*y3[i]))-y3[i]**2).real
#    f_out.write(str(x)+" "+str(y)+'\n')
#    f_out2.write(str(x)+" "+str(y)+'\n')
##f_out.close()    # () at the end is necessary to close the file
##f_out2.close()    # () at the end is necessary to close the file
#y4 = numpy.linspace(-0.7,0,n)
##y4 = numpy.linspace(-0.062,0,n)      
#for i in range(len(y4)):    
#    y = y4[i]
#    x = cmath.sqrt((1.0/4)+(y4[i]/math.tan(c1*y4[i]))-y4[i]**2).real
#    f_out.write(str(x)+" "+str(y)+'\n')
#    f_out2.write(str(x)+" "+str(y)+'\n')
f_out.close()    # () at the end is necessary to close the file
f_out1.close()    # () at the end is necessary to close the file
f_out2.close()    # () at the end is necessary to close the file