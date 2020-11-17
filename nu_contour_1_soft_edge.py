# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 09:14:32 2018

@author: yadav
"""
"""
Created on Fri Nov 16 09:50:09 2018

@author: yadav
"""

# code to generate contour nu, nu1 and nu2 for J(s)=c1s+c0-log((s-1/2)/(s+1/2))
#V(x)=x^2/2 ; so c1=1,c0=1/2
#Due to given from of locus equation we first find max and min values of y(here for x=0)
#then we vary y such that we get counterclockwise contour
#stepsize for each quarter of plot can be chosen
import math
import cmath
import numpy
import matplotlib    #pylab is submodule in matplotlib

c1 = 2.0
c1_short = 2.0
c0 = 1.0
c0_short = 1.0
iteration = 1.0
#c1 = 50.0    # to check hard edge behaviour as c1->inf so let c1=50,c0=25, so V(x)=x^2/100
#c0 = 25.0   
n = 10000    # number of points in each quadrant
f_out=file("input_output_files/rho_2.0/contour/nu_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_"+str(n)+"points_iter"+str(iteration)+".txt","w")
f_out1=file("input_output_files/rho_2.0/contour/nu1_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_"+str(n)+"points_iter"+str(iteration)+".txt","w")
f_out2=file("input_output_files/rho_2.0/contour/nu2_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_"+str(n)+"points_iter"+str(iteration)+".txt","w")
y1 = numpy.linspace(0,1.7,n)    #linspace creates array of x values, y_max=1.05 for c1=1,c0=1/2  
#y1 = numpy.linspace(0,0.062,n)    #linspace creates array of x values, y_max=0.062 for c1=50,c0=25 
#x = numpy.empty(len(y1),float)
#y = numpy.empty(len(y1),float)
for i in range(len(y1)):    # function len() on array gives no. of rows of array
    y = y1[i]
    x = cmath.sqrt((1.0/4)+(y1[i]/math.tan(c1*y1[i]))-y1[i]**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative
    f_out.write(str(x)+" "+str(y)+'\n')
    f_out1.write(str(x)+" "+str(y)+'\n')
#f_out.close()    # () at the end is necessary to close the file 
#f_out1.close()    # () at the end is necessary to close the file
y2 = numpy.linspace(1.7,0,n)
 
#y2 = numpy.linspace(0.062,0,n)     
for i in range(len(y2)):    
    y = y2[i]
    x = (-1)*cmath.sqrt((1.0/4)+(y2[i]/math.tan(c1*y2[i]))-y2[i]**2).real
    f_out.write(str(x)+" "+str(y)+'\n')
    f_out1.write(str(x)+" "+str(y)+'\n')
#f_out.close()    # () at the end is necessary to close the file 
#f_out1.close()    # () at the end is necessary to close the file
y3 = numpy.linspace(0,-1.7,n)
#y3 = numpy.linspace(0,-0.062,n)      
for i in range(len(y3)):    
    y = y3[i]
    x = (-1)*cmath.sqrt((1.0/4)+(y3[i]/math.tan(c1*y3[i]))-y3[i]**2).real
    f_out.write(str(x)+" "+str(y)+'\n')
    f_out2.write(str(x)+" "+str(y)+'\n')
#f_out.close()    # () at the end is necessary to close the file
#f_out2.close()    # () at the end is necessary to close the file
y4 = numpy.linspace(-1.7,0,n)
#y4 = numpy.linspace(-0.062,0,n)      
for i in range(len(y4)):    
    y = y4[i]
    x = cmath.sqrt((1.0/4)+(y4[i]/math.tan(c1*y4[i]))-y4[i]**2).real
    f_out.write(str(x)+" "+str(y)+'\n')
    f_out2.write(str(x)+" "+str(y)+'\n')
f_out.close()    # () at the end is necessary to close the file
f_out1.close()    # () at the end is necessary to close the file
f_out2.close()    # () at the end is necessary to close the file
#print len(y1)+len(y2)+len(y3)+len(y4)
data = numpy.loadtxt("input_output_files/rho_2.0/contour/nu_contour_c1="+str(c1_short)+"_c0="+str(c0_short)+"_"+str(n)+"points_iter"+str(iteration)+".txt",float)
x = data[:,0]
y = data[:,1]
matplotlib.pylab.plot(x,y)
matplotlib.pylab.ylim(-1.15,1.15)    
