clc; close all; clear all;


% Set system parameters 
l_1 = 1;            % length of upper rod (rod 1)
l_2 = 1;            % length of lower rod (rod 2)
m_1 = 5;            % mass of upper rod (rod 1)
m_2 = 5;            % mass of lower rod (rod 2)
g = 9.81;           % gravity 

% calculate mass moment of inertias 
I_1 = 1/12*m_1*l_1^2;       % mass moment of inertia of rod 1
I_2 = 1/12*m_2*l_2^2;       % mass moment of inertia of rod 2

% set initial condition of rod 1 
theta_1 = pi/4;
x_1 = l_1/2*sin(theta_1);
y_1 = -l_1/2*cos(theta_1);
vx_1 = 0;
vy_1 = 0;
theta_dot_1 = 0;

% set initial conditions of rod 2
theta_2 = pi/4;
x_2 = l_2/2*sin(theta_2);
y_2 = -l_2/2*cos(theta_2);
vx_2 = 0;
vy_2 = 0;
theta_dot_2 = 0;

% define Mass metrix
M(6:6) = 0;     % Initilize matrix
M(1,1) = m_1; M(2,2) = m_1; M(3,3) = I_1;
M(4,4) = m_2; M(5,5) = m_2; M(6,6) = I_2;

% define force vector
h = [0;-m_1*g;0;0;-m_2*g;0];

% set parameters for neumerical analysis
n = 1000000;        % number of points for the analysis
time_limit = 20;    % simulation time limit

%set time array
t(1:n+1)=0;
for i=1:n+1
    t(i)= time_limit/n*(i-1);
end
% set the stepsize
dt= t(2)-t(1);

% initilize the velocity and position vectors
v(1:6,n+1) = 0;
r(1:6,n+1) = 0;

% add initial condition on the velocity and position vectors
v(1:6,1) = [vx_1;vy_1;theta_dot_1;vx_2;vy_2;theta_dot_2];
r(1:6,1) = [x_1;y_1;theta_1;x_2;y_2;theta_2];

% initilize jacobiam and derivative of the jacobian
D(1:4,1:6) = 0;
D_dot(1:4,1:6) = 0;

% Initilize the loop for numerical intigeration

for i=2:n+1
    



end








