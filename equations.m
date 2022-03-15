%%%%%% adapted from equations_falling_block.m %%%
clc
clear all 
close all

syms x y z real
syms vx vy vz real
syms ax ay az real
syms phi theta psi real %euler angles
syms phidot thetadot psidot real
syms phiddot thetaddot psiddot real 
syms m g Ixx Iyy Izz real 
syms K l b Ax Ay Az real
syms omega1 omega2 omega3 omega4 real

i = [1 0 0]';
j = [0 1 0]';
k = [0 0 1]';

%%%%%%%%%% Rotation Matrices %%%%%%%%%%
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
   
%%%%% angular velocity and energy %%%%
om_b = phidot*i +  R_x'*(thetadot*j) + R_x'*R_y'*(psidot*k); % angular velocities in body frame
I = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];%body frame inertia[Diagonal Matrix]
v = [vx; vy; vz];
T = 0.5*m*(v')*v + 0.5*om_b'*I*om_b;
V = m*g*z;
L = T-V;
% disp('copy paste energy in the code');
disp(['KE(i) = ',char(T),';']);
disp(['PE(i) = ',char(V),';']);
disp(['TE(i) = KE(i)+PE(i);']);
disp(' ');

%Derive equations of motion using Euler-Lagrange method
q = [x y z  phi theta psi ];
qdot = [vx vy vz phidot thetadot psidot ];
qddot = [ax ay az phiddot thetaddot psiddot];


%%%%%% external forces and torques %%%%%%%%%
Thrust = [0; 0; K*(omega1^2+omega2^2+omega3^2+omega4^2)];
Drag = [Ax*vx; Ay*vy; Az*vz];
F_ext = R*Thrust-Drag;
tau_phi = K*l*(omega4^2 - omega2^2);
tau_theta = K*l*(omega3^2 - omega1^2);
tau_psi = b*(omega1^2-omega2^2+omega3^2-omega4^2);
tau_ext = [tau_phi; tau_theta; tau_psi];

T_ext = [F_ext; tau_ext];



for ii=1:6
    dLdqdot(ii) = diff(L,qdot(ii));
    ddt_dLdqdot(ii) = diff(dLdqdot(ii),q(1))*qdot(1) + diff(dLdqdot(ii),qdot(1))*qddot(1)+...
                      diff(dLdqdot(ii),q(2))*qdot(2) + diff(dLdqdot(ii),qdot(2))*qddot(2)+...
                      diff(dLdqdot(ii),q(3))*qdot(3) + diff(dLdqdot(ii),qdot(3))*qddot(3)+...
                      diff(dLdqdot(ii),q(4))*qdot(4) + diff(dLdqdot(ii),qdot(4))*qddot(4)+...
                      diff(dLdqdot(ii),q(5))*qdot(5) + diff(dLdqdot(ii),qdot(5))*qddot(5)+...
                      diff(dLdqdot(ii),q(6))*qdot(6) + diff(dLdqdot(ii),qdot(6))*qddot(6);
    dLdq(ii) = diff(L,q(ii));

    EOM(ii) = ddt_dLdqdot(ii) - dLdq(ii) - T_ext(ii);
end


%%%%%%% post process equations %%%%%%%
A = jacobian(EOM,qddot);
for ii=1:6
    B(ii,1) = -subs(EOM(ii),qddot,[0 0 0 0 0 0]);
end

disp('copy paste in MATLAB');
disp(' ');
for ii=1:6
    for jj=1:6
        string = [ 'A(',num2str(ii),',',num2str(jj),')='];
        disp([string, char(simplify(A(ii,jj))), ';']);
    end
end

disp(' ');
for ii=1:6
    string = [ 'B(',num2str(ii),',1)='];
    disp([string, char(simplify(B(ii,1))), ';']);
end

disp(' ');
disp('X = A\B;');
disp(' ');