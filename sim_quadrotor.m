function sim_quadcopter
%%%%%% adapted from sim_falling_block.m %%%
clc
clear all
close all

%Intializing mass, inertia, gravity, drag coefficients
params.m = 0.468;
params.g = 9.81;
params.Ixx = 4.856*1e-3;
params.Iyy = 4.856*1e-3;
params.Izz = 8.801*1e-3;
params.l = 0.225;
params.K = 2.980*1e-6;
params.b = 1.14*1e-7;
params.Ax = 0.25*0;
params.Ay = 0.25*0;
params.Az = 0.25*0;

%%%%%%%%%% play with speeds of each motor to get the quadcopter to roll,
%%%%%%%%%% yaw, pitch, and hover  %%%%%%%%%%
%%%%% Tz = K*(omega1^2+omega2^2+omega3^2+omega4^2); balances with mg
%%%%% tau_x = tau_phi = K*L*(omega4^2 - omega2^2);
%%%%% tau_y = tau_theta = K*L*(omega3^2 - omega1^2);
%%%%% tau_z = tau_psi = b*(omega1^2-omega2^2 + omega3^2-omega4^2);

omega = 1.075; %1.075
speed = omega*sqrt(1/params.K);
dspeed = 0.05*speed;
%params.omega1 = speed;
%params.omega2 = speed;
%params.omega3 = speed;
%params.omega4 = speed;

%%% to turn about y-axis
 %dspeed = 0.05*speed;
 %params.omega1 = speed;
 %params.omega2 = speed;
 %params.omega3 = speed-dspeed;
 %params.omega4 = speed;

%%% to spin about the x axis
% dspeed = 0.05*speed;
% params.omega1 = speed;
% params.omega2 = speed-dspeed;
% params.omega3 = speed;
% params.omega4 = speed;

%%% to hover and spin about z axis
params.omega1 = speed-0.5*dspeed;
params.omega2 = speed+0.5*dspeed;
params.omega3 = speed-0.5*dspeed;
params.omega4 = speed+0.5*dspeed;


%position and velocity (6 degrees of freedom)
x0 = 0;
y0 = 0;
z0 = 0;
vx0 = 0;
vy0 = 0;
vz0 = 0;
	
% initialize euler angles and rates
phi0 = 0;
theta0 = 0;
psi0 = 0;
phidot0 = 0;
thetadot0 = 0; 
psidot0 = 0;

%Integrate with ODE
tspan = linspace(0,1); %time span
Z0=[x0 y0 z0 phi0 theta0 psi0 vx0 vy0 vz0 phidot0 thetadot0 psidot0]';
options=odeset('abstol',1e-12,'reltol',1e-12);

[T Z] = ode45(@eom,tspan,Z0,options,params); % This ODE trying to take the eqaution of motion


figure(1)
plot(T,Z);
legend('x','y','z','$\phi$','$\theta$','$\psi$','$\dot{x}$','$\dot{y}$','$\dot{z}$','$\dot{\phi}$','$\dot{\theta}$','$\dot{\psi}$','Interpreter','latex','Fontsize',12)

figure(2)
positions = Z(:,1:3);
angles = Z(:,4:6);
animate(positions,angles,params.l)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Zdot = eom(t,Z,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = params.m;
Ixx = params.Ixx;
Iyy = params.Iyy;
Izz = params.Izz;
g = params.g;
l = params.l;
K = params.K;
b = params.b;
Ax = params.Ax;
Ay = params.Ay;
Az = params.Az;

omega1 = params.omega1;
omega2 = params.omega2;
omega3 = params.omega3;
omega4 = params.omega4;

x = Z(1); y = Z(2); z = Z(3);
phi = Z(4); theta = Z(5); psi = Z(6);
vx = Z(7); vy = Z(8); vz = Z(9);
phidot = Z(10); thetadot = Z(11); psidot = Z(12);


%%%%%%% copy pasted from running equations code in matlab %%%%%%%%%
A(1,1)=m;
A(1,2)=0;
A(1,3)=0;
A(1,4)=0;
A(1,5)=0;
A(1,6)=0;
A(2,1)=0;
A(2,2)=m;
A(2,3)=0;
A(2,4)=0;
A(2,5)=0;
A(2,6)=0;
A(3,1)=0;
A(3,2)=0;
A(3,3)=m;
A(3,4)=0;
A(3,5)=0;
A(3,6)=0;
A(4,1)=0;
A(4,2)=0;
A(4,3)=0;
A(4,4)=Ixx;
A(4,5)=0;
A(4,6)=-Ixx*sin(theta);
A(5,1)=0;
A(5,2)=0;
A(5,3)=0;
A(5,4)=0;
A(5,5)=Iyy - Iyy*sin(phi)^2 + Izz*sin(phi)^2;
A(5,6)=cos(phi)*cos(theta)*sin(phi)*(Iyy - Izz);
A(6,1)=0;
A(6,2)=0;
A(6,3)=0;
A(6,4)=-Ixx*sin(theta);
A(6,5)=cos(phi)*cos(theta)*sin(phi)*(Iyy - Izz);
A(6,6)=Ixx*sin(theta)^2 + Izz*cos(phi)^2*cos(theta)^2 + Iyy*cos(theta)^2*sin(phi)^2;
 
B(1,1)=K*(sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(omega1^2 + omega2^2 + omega3^2 + omega4^2) - Ax*vx;
B(2,1)=- Ay*vy - K*(cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(omega1^2 + omega2^2 + omega3^2 + omega4^2);
B(3,1)=K*cos(phi)*cos(theta)*(omega1^2 + omega2^2 + omega3^2 + omega4^2) - g*m - Az*vz;
B(4,1)=(Izz*thetadot^2*sin(2*phi))/2 - (Iyy*thetadot^2*sin(2*phi))/2 - K*l*omega2^2 + K*l*omega4^2 + Ixx*psidot*thetadot*cos(theta) - Iyy*psidot*thetadot*cos(theta) + Izz*psidot*thetadot*cos(theta) + Iyy*psidot^2*cos(phi)*cos(theta)^2*sin(phi) - Izz*psidot^2*cos(phi)*cos(theta)^2*sin(phi) + 2*Iyy*psidot*thetadot*cos(phi)^2*cos(theta) - 2*Izz*psidot*thetadot*cos(phi)^2*cos(theta);
B(5,1)=(Ixx*psidot^2*sin(2*theta))/2 - K*l*omega1^2 + K*l*omega3^2 - Ixx*phidot*psidot*cos(theta) + Iyy*phidot*thetadot*sin(2*phi) - Izz*phidot*thetadot*sin(2*phi) - Izz*psidot^2*cos(phi)^2*cos(theta)*sin(theta) - Iyy*psidot^2*cos(theta)*sin(phi)^2*sin(theta) - Iyy*phidot*psidot*cos(phi)^2*cos(theta) + Izz*phidot*psidot*cos(phi)^2*cos(theta) + Iyy*phidot*psidot*cos(theta)*sin(phi)^2 - Izz*phidot*psidot*cos(theta)*sin(phi)^2;
B(6,1)=b*omega1^2 - b*omega2^2 + b*omega3^2 - b*omega4^2 + Ixx*phidot*thetadot*cos(theta) + Iyy*phidot*thetadot*cos(theta) - Izz*phidot*thetadot*cos(theta) - Ixx*psidot*thetadot*sin(2*theta) + Iyy*psidot*thetadot*sin(2*theta) + Iyy*thetadot^2*cos(phi)*sin(phi)*sin(theta) - Izz*thetadot^2*cos(phi)*sin(phi)*sin(theta) - 2*Iyy*phidot*thetadot*cos(phi)^2*cos(theta) + 2*Izz*phidot*thetadot*cos(phi)^2*cos(theta) - 2*Iyy*phidot*psidot*cos(phi)*cos(theta)^2*sin(phi) + 2*Izz*phidot*psidot*cos(phi)*cos(theta)^2*sin(phi) - 2*Iyy*psidot*thetadot*cos(phi)^2*cos(theta)*sin(theta) + 2*Izz*psidot*thetadot*cos(phi)^2*cos(theta)*sin(theta);
 
X = A\B;

Zdot = [vx vy vz phidot thetadot psidot X(1) X(2) X(3) X(4) X(5) X(6)]'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function animate(positions,angles,l)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axle_x = [-l/2 0 0;
           l/2 0 0];
axle_y = [0 -l/2 0;
          0  l/2 0];
      
dmax = max([max(positions),l]);

r = 0.1*l; %radius of propellers
ang = linspace(0,2*pi);
x_circle = r*cos(ang);
y_circle = r*sin(ang);
z_circle = zeros(1,length(ang));
propeller = [x_circle',y_circle',z_circle'];


[p1,q1] = size(propeller);
[p2,q2] = size(axle_x);
[mm,nn] = size(angles);

for ii=1:mm
    x = positions(ii,1);
    y = positions(ii,2);
    z = positions(ii,3);
    phi = angles(ii,1); 
    theta = angles(ii,2);
    psi = angles(ii,3);
    R = get_rotation(phi,theta,psi);
    
    for i=1:p2
        r_body = axle_x(i,:)';
        r_world = R*r_body;
        new_axle_x(i,:) = r_world';
    end
    new_axle_x = [x y z] +new_axle_x;
    
    for i=1:p2
        r_body = axle_y(i,:)';
        r_world = R*r_body;
        new_axle_y(i,:) = r_world';
    end
    new_axle_y = [x y z] +new_axle_y;
    
    for i=1:p1
        r_body = propeller(i,:)';
        r_world = R*r_body;
        new_propeller(i,:) = r_world'; 
    end
    new_propeller1 = new_axle_x(1,:) + new_propeller;
    new_propeller3 = new_axle_x(2,:) + new_propeller;
    new_propeller2 = new_axle_y(1,:) + new_propeller;
    new_propeller4 = new_axle_y(2,:) + new_propeller;

     line(new_axle_x(:,1),new_axle_x(:,2),new_axle_x(:,3),'Linewidth',2); hold on;
     line(new_axle_y(:,1),new_axle_y(:,2),new_axle_y(:,3),'Linewidth',2);
    patch(new_propeller1(:,1),new_propeller1(:,2),new_propeller1(:,3),'r');
    patch(new_propeller2(:,1),new_propeller2(:,2),new_propeller2(:,3),'g');
    patch(new_propeller3(:,1),new_propeller3(:,2),new_propeller3(:,3),'b');
    patch(new_propeller4(:,1),new_propeller4(:,2),new_propeller4(:,3),'c');
    axis(1.2*[-dmax dmax -dmax dmax -dmax dmax]);
    xlabel('x'); ylabel('y'); zlabel('z');
    %view(0,90)
    view(3)
    pause(0.01)
    if (ii~=mm)
        clf
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = get_rotation(phi,theta,psi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%% uses 3-2-1 euler angles
R_x = [1    0       0; ...
       0  cos(phi) -sin(phi); ...
       0  sin(phi) cos(phi)];
   
R_y = [cos(theta)  0   sin(theta); ...
       0           1         0; ...
      -sin(theta) 0   cos(theta)]; 
   
R_z = [cos(psi) -sin(psi)  0; ...
       sin(psi)  cos(psi)  0; ...
       0           0       1];  
R = R_z*R_y*R_x;
