# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 17:40:59 2019

@author: yadav
"""

# code to generate contour nu, nu1 and nu2 for J(s)=(c1s+c0)((s+1)/s)**(1/theta)
import math
import cmath
import numpy

#rho = 2.0    # V(x)=rho*x
#c = theta/rho    # eqn.(4.21) in C.R.


#def contour_CR_prob_soft_edge(c1,c0):
c1=2.0
c0=1.0
n = 100000    # number of points in each quadrant
    
y1 = numpy.linspace(0.7/n,0.7,n)    #linspace creates array of x values, y_max=1.05 for c1=1,c0=1/2  
f_out=file("input_output_files/rho_2.0/contour/nu_contour_function_soft_edge_10000points.txt","w")
#    f_out=open("input_output_files//contour//nu_contour_function_soft_edge_18000points.txt","w")    # open() can also be used to create file    
y1 = numpy.linspace(0,0.7,n)    #linspace creates array of x values, y_max=1.05 for c1=1,c0=1/2  

for i in range(len(y1)):    # function len() on array gives no. of rows of array
    y = y1[i]
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0:
        x = cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real    # we take cmath.sqrt() because sometimes argument of sqrt() is negative
        f_out.write(str(x)+" "+str(y)+'\n')

y2 = numpy.linspace(0.7,0.7/n,n)

for i in range(len(y2)):    
    y = y2[i]
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0:
        x = (-1)*cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real
        f_out.write(str(x)+" "+str(y)+'\n')
        
y3 = numpy.linspace(-0.7/n,-0.7,n)
    
for i in range(len(y3)):    
    y = y3[i]
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0:
        x = (-1)*cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real
        f_out.write(str(x)+" "+str(y)+'\n')
        
y4 = numpy.linspace(-0.7,-0.7/n,n)
    
for i in range(len(y4)):    
    y = y4[i]
    if ((1.0/4)+(y/math.tan(c1*y))-y**2)>=0 and i<=9000:
        x = cmath.sqrt((1.0/4)+(y/math.tan(c1*y))-y**2).real
        f_out.write(str(x)+" "+str(y)+'\n')
        
f_out.close()    # () at the end is necessary to close the file

#    data = numpy.loadtxt("input_output_files/rho_2.0/contour/nu_contour_function_soft_edge_10000points.txt",float)      
#    contour = data[:,0:2] # The plain [:] operator slices from beginning to end-1, Eg. a[:,0:n] gives all the rows from 0th column to n-1 column    
    
#    return contour                
