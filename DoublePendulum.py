import numpy as np
import sympy as smp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.animation import PillowWriter

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

#define the equation of the Potential, Keinetic Energy and Lagrangian
T= 1/2 *m1* (smp.diff(x1,t)**2 + smp.diff(y1,t)**2)+\
1/2 *m2* (smp.diff(x2,t)**2 + smp.diff(y2,t)**2)
V = (m1*g*y1) + (m2*g*y2) + 1/2 * k * r1** 2
L= T-V

#calculate Lagrangian’s equation and solve for equation
#of motion of the system
LE1 = smp.diff(L, theta1) - smp.diff(smp.diff(L, theta1_d), t)
LE1 = LE1.simplify()

LE2 = smp.diff(L, r1) - smp.diff(smp.diff(L, r1_d), t)
LE2 = LE2.simplify()

LE3 = smp.diff(L, theta2) - smp.diff(smp.diff(L, theta2_d), t)
LE3 = LE3.simplify()

sols = smp.solve([LE1,LE2,LE3],(theta1_dd,r1_dd,theta2_dd),\
simplify=True, rational = True)

#reduction of 2nd order to 1st order ODE
dw1dt_f = smp.lambdify((k,g,m1,l1,m2,l2,theta1,theta1_d,r1,r1_d,theta2,theta2_d),\
sols[theta1_dd],modules=['numpy'])
dtheta1dt_f = smp.lambdify(theta1_d,theta1_d,modules=['numpy'])

dv1dt_f = smp.lambdify((k,g,m1,l1,m2,l2,theta1,theta1_d,r1,r1_d,theta2,theta2_d),\
sols[r1_dd],modules=['numpy'])
dr1dt_f = smp.lambdify(r1_d,r1_d,modules=['numpy'])

dw2dt_f = smp.lambdify((k,g,m1,l1,m2,l2,theta1,theta1_d,r1,r1_d,theta2,theta2_d),\
sols[theta2_dd],modules=['numpy'])
dtheta2dt_f = smp.lambdify(theta2_d,theta2_d,modules=['numpy'])

#create a function to calculate the total energy
E = T+V
E = E.simplify()
Energy= smp.lambdify((k,g,m1,l1,m2,l2,theta1,theta1_d,r1,r1_d,theta2,theta2_d)\
,E,modules=['numpy'])
            
#assemble all of the equations in one function

def system(variables, t):
    theta1,w1,r1,v1,theta2,w2 = variables
    return [
    dtheta1dt_f(w1),
    dw1dt_f(k,g,m1,l1,m2,l2,theta1,w1,r1,v1,theta2,w2),
    dr1dt_f(v1),
    dv1dt_f(k,g,m1,l1,m2,l2,theta1,w1,r1,v1,theta2,w2),
    dtheta2dt_f(w2),
    dw2dt_f(k,g,m1,l1,m2,l2,theta1,w1,r1,v1,theta2,w2)
    ]

# define the System parameters
g=9.81
k=100
m1=2
m2=2
l1=2
l2=2
# define Initial Conditions
theta1_0 = 5*np.pi/4
theta1_d_0 = 0
r1_0 = 0
r1_d_0 = 0
theta2_0 = 5*np.pi/4
theta2_d_0 = 0

#Numeric Solution parameters
h= 0.001 # step size
simulation_time = 10 # simulation time in seconds
grids = int(simulation_time/h)
grid_points = np.linspace(0,simulation_time,grids)

y0 = np.array([theta1_0,theta1_d_0,r1_0,r1_d_0,theta2_0,theta2_d_0])
ans = np.zeros((len(grid_points),6)) #
ans[0]=y0
t=0

# Numeric Solvers
for i in range(1,len(grid_points)):

# Forward Euler method
   # y0=list(np.array(y0))+(np.array(system(ans[i-1],t))*h)

   # ans[i]= y0
   # t=t+h


# Predictor Corrector Method
   # predict=list(np.array(y0))+(np.array(system(ans[i-1],t))*h)
   # correct = list(np.array(y0))+(np.array(system(predict,t))*h)
   # y0= (predict+correct)/2
   # ans[i]= y0
   # t=t+h



# Runge-Kutta method
   k1 = system(y0,t)
   k2=system(list(np.array(y0)+(np.array(k1)*h/2)),t+h/2)
   k3=system(list(np.array(y0)+(np.array(k2)*h/2)),t+h/2)
   k4=system(list(np.array(y0)+(np.array(k3)*h)),t+h)
   y0=list(np.array(y0)+((np.array(k1)+np.array(k2)*2+np.array(k3)*2+np.array(k4)))*h/6)
   ans[i]= y0
   t=t+h

# plotting
def get_coord(l1,l2,theta1 , r1, theta2):
   return ((l1+r1)*np.cos(theta1),-(l1+r1)*np.sin(theta1),(l1+r1)*np.cos(theta1)+l2*np.cos(theta2),-(l1+r1)*np.sin(theta1)-l2*np.sin(theta2))


x1, y1, x2, y2 = get_coord(l1,l2,ans.T[0], ans.T[2], ans.T[4])

trajectory_x1=[]
trajectory_y1=[]
trajectory_x2=[]
trajectory_y2=[]

def animate(i):
   i= i * int(grids/(simulation_time*50))
   ln1.set_data([0, x1[i]], [0, y1[i]])
   ln2.set_data([x1[i], x2[i]], [y1[i], y2[i]])
   trajectory_x1.append(x1[i])
   trajectory_y1.append(y1[i])
   trajectory_x2.append(x2[i])
   trajectory_y2.append(y2[i])
   trajectory1.set_data(trajectory_x1,trajectory_y1)
   trajectory2.set_data(trajectory_x2,trajectory_y2)
   return ln1,ln2,trajectory1,trajectory2

fig, ax = plt.subplots(1,1, figsize=(8,8))
ax.grid()
ln1, = plt.plot([], [], 'ro--', lw=3, markersize=1)
ln2, = plt.plot([], [], 'ro-', lw=3, markersize=8)
trajectory1, = ax.plot([], [],'-',color='green',alpha=0.5)
trajectory2, = ax.plot([], [],'-',color='blue',alpha=0.5)
ax.set_ylim(-8, 8)
ax.set_xlim(-8,8)
ani = animation.FuncAnimation(fig, animate, frames=int(simulation_time*50), interval =50 )

ani.save('spring_pendulum.gif',writer='pillow',fps=50)
fig.savefig('final_frame.png')

Eng=[]
for i in range(0,len(grid_points)):
    Eng.append(Energy(k,g,m1,l1,m2,l2,ans.T[0,i],ans.T[1,i],ans.T[2,i],ans.T[3,i],ans.T[4,i],ans.T[5,i]))

figure2 = plt.figure(2)
figure2 =plt.plot(grid_points,Eng)
figure2=plt.show()

delta_TE=np.zeros((grids, 1))
for i in range(0,len(Eng)-1) :
    delta_TE[i,:]=np.array(Eng[0]-Eng[i])

figure3 = plt.figure(3)
figure3 = plt.plot(grid_points,delta_TE)
figure3 = plt.show()

print('maximum error =')
print(float(max(abs(delta_TE))))