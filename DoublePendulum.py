import numpy as np
import sympy as smp

# define variables constant variables
t,m1,l1,m2,l2,g,k = smp.symbols('t m_1 l_1 m_2 l_2 g k')
# define changing variables (degree of freedom)
theta1,theta2,r1 = smp.symbols('theta_1, theta_2, r_1', cls = smp.Function)
#make them function of time
theta1 = theta1(t)
r1= r1(t)
theta2= theta2(t)
#define first derivative
theta1_d = smp.diff(theta1,t)
r1_d = smp.diff(r1,t)
theta2_d = smp.diff(theta2,t)
# define second derivative
theta1_dd = smp.diff(theta1_d,t)
r1_dd = smp.diff(r1_d,t)
theta2_dd = smp.diff(theta2_d,t)
# define coordniate variables
x1,y1,x2,y2= smp.symbols('x_1,y_1,x_2,y_2', cls = smp.Function)
# make them function of respective variables
x1 = x1(theta1)
y1 = y1(theta1)
x2 = x2(theta1,r1,theta2)
y2 = y2(theta1,r1,theta2)
# define coordinates of the mass
x1= (l1+r1)*smp.cos(theta1)
y1= -(l1+r1)*smp.sin(theta1)
x2= x1 + l2*smp.cos(theta2)
y2= y1 - l2*smp.sin(theta2)
