clear all; clc;

roll = deg2rad(12); pitch = deg2rad(24); yaw = deg2rad(36);
rpy = [roll, pitch, yaw];
R = rpyToRot(rpy);
[rpy2, rpy3] = rotToRpy(R);

assert( isEqual(R.'*R, eye(3) ))
assert( isEqual(rpy3.' , rpy))

eul_rate_B = (1/cos(pitch)) * [0 ,          sin(roll),              cos(roll);
                               0 ,          cos(roll)*cos(pitch),   -sin(roll)*cos(pitch);
                              cos(pitch),   sin(roll)*sin(pitch),    cos(roll)*sin(pitch)];
                          
omega_body = rand(3,1);
eul_rate = eul_rate_B * omega_body;

%% Test Euler-Rate
dt = 1e-4;

Rnew = R * expm(skew(omega_body)*dt); % R * exp(omegax dt)
[rpy2new, rpy3new] = rotToRpy(Rnew);

droll = rpy3new(1) - roll;
dpitch = rpy3new(2) - pitch;
dyaw = rpy3new(3) - yaw;

eul_rate_num = [dyaw, dpitch,droll ].'/dt;

tol = 1e-4;  
diffVec = eul_rate_num - eul_rate;

assert(norm(diffVec) < tol, ...
    'eul_rate_B test failed: difference is too large!')
fprintf('Test passed: eul_rate_B is consistent with time-stepping. ||Diff|| = %f\n', ...
    norm(diffVec));

%% 
syms phi theta psi % roll, pitch, yaw

Rsym = rpyToRot([phi, theta, psi]);
dR_dphi = diff(Rsym, phi);
dR_dtheta = diff(Rsym, theta);
dR_dyaw = diff(Rsym, psi);

yaw_dot = eul_rate(1);
pitch_dot = eul_rate(2);
roll_dot = eul_rate(3);

expr = R*skew(omega_body); % R omega x
expr2 = double(subs(dR_dphi, [phi, theta, psi],[roll, pitch, yaw])*roll_dot+...
              subs(dR_dtheta, [phi, theta, psi],[roll, pitch, yaw])*pitch_dot+...
              subs(dR_dyaw, [phi, theta, psi],[roll, pitch, yaw])*yaw_dot);

assert( isEqual(expr, expr2))



