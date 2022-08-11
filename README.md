# UGRP

#Program overview
This program can simulate ray dynamics in a micro-cavity consists of 3 partial ellipses.
There are 5 output data files: boundary.dat, ray.dat, psos.dat, spd.dat and ffield.dat which means that with this program, we can draw **ray trajectory**, **poincare surface of section, survival probability distribution and far-field pattern.**
Language: Fortran90
Environment: Ubuntu LTS 20.04
Plotting: GNUPlot

#Microcavity
A microcavity consists of three partial ellipses.
A half ellipse and two quarter ellipses. 
Half ellipse: (x/g1)^2 + y^2 = 1
quarter ellipse1: ((x-(1-g2))/(2-g2))^2 + (y/g3)^2 = 1
quarter ellipse2: ((x-(1-g2))/g2)^2 + (y/g3)^2 = 1

#Hyper parameter description
There are **11 hyper parameters** including g1, g2 and g3.
**g1, g2 and g3** are kind of deformation parameters related to the shape of each ellipse.
**s_num and p_num** are parameters that determine initial position and angle of the ray. 
s_num is the number that divide 360 degree ex) s_num = 3 -> initial position: 120, 240, 360 in degree
p_num is the number that divide value of sin (-1 to 1) ex) p_num = 4 -> sin(x) = -0.5, 0, 0.5 -> -30, 0, 30 in degree 
**num** is the number of iterations of each ray.
**n_i and n_t** are refraction index of microcavity and its surrounding respectively
**row and col** are the number of grid of SPD.
**l** is the observation distance of far-field pattern.

#Problems 1
The program takes too many times to execute. 
Innately, it requires enormous calculation (In my scale). Ex) numerically calculate the arc length every ray reflection.
And writing data file may take long time.
Regardless of the reasons above, I think my algorithm makes the program much more slower. 
There might be several reasons but I don't know what actually they are.
The only thing that come to my mind is a lot of if statement.
Please let me know if you find some slow-down factors. I'll deeply appreciate that.

#Problems 2
When it comes to calculatation of 's' which is x-axis of phase space, the calculation accuracy is not high enough.
I adopt trapezoidal rule to numerically caluculate the arc length integral.

#Probelms 3
PSOS are imperfect
The section of quarter ellipse2 has more data than other sections. And some 's' values are slightly larger than 1 because of 's' calculation error 


