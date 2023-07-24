function [robot, params] = build_quadruped_model(mode)
% coordinate convention
% x+ moving direction
% z+ vetical upward
% y+ perp to x and z and agree with right hand rule
% hip link and knee link are on the negative direction of the corresponding
% link-attached frame
% all rotational inertias are w.r.t. CoM, which is required by mcI 

% in 2D quadruped case, mass, CoM and rotational inertia of the body are
% updated to incoporate the four abads

%% get mini-cheetah physical parameters (from Cheetach Software)
params = get_miniCheetah_params();
params.g = 9.81;

CoM_body = params.CoM_body;
CoM_hip = params.CoM_hip;
CoM_knee = params.CoM_knee; 
CoM_abad = params.CoM_abad;

mass_body = params.mass_body;
mass_hip = params.mass_hip;
mass_knee = params.mass_knee;
mass_abad = params.mass_abad;
params.mass_robot = mass_body + 4*(mass_hip + mass_knee + mass_abad);

% link rotational inertia w.r.t. its CoM expressed in local frame
I_body = params.I_body;     % trunk rotational inertia
I_hip = params.I_hip;       % hip rotational inertia
I_knee = params.I_knee;     % knee rotational inertia
I_abad = params.I_abad;     % abad rotational inertia

bodyLength = params.bodyLength;
bodyWidth = params.bodyWidth;
hiplinkLength = params.hiplinkLength;
hip_x = [bodyLength -bodyLength]/2;
params.hip_x = hip_x;       % used in forward kinematics

loc_abad = [bodyLength, bodyWidth, 0;
            bodyLength, -bodyWidth, 0;
            -bodyLength, bodyWidth, 0;
            -bodyLength, -bodyWidth, 0]'/2; % locations of four abads

%% generate 2D-quadruped model
if strcmp(mode,'2D')        % double the mass inertia in 2D quadruped
    NLEGS = 2;              % two legs in planar quadruped (NLEGS = 4 for 3D quadruped)
    mass_hip = mass_hip * 2;
    mass_knee = mass_knee*2;
    mass_body = mass_body + 4*mass_abad;    % 4 abads in body
    
    % update body CoM (considering 4 abads into body)
    weighted_loc = 0;
    for i = 1:4
        weighted_loc = weighted_loc + mass_abad*(loc_abad(:,i)+CoM_abad);
    end
    CoM_body = weighted_loc/mass_body;
        
    I_hip = I_hip*2;        % two legs in one
    I_knee = I_knee*2;      % two legs in one
    
    % update body inertia
    S = skew(CoM_body); % skew(v) defined in spatial v2
    I_body = I_body + mass_body*(S*S');     % body inertia w.r.t. new CoM
    for i = 1:4
        S = skew(loc_abad(:,i) + CoM_abad - CoM_body);
        I_body = I_body + I_abad + mass_abad*(S*S');    % body inertia considering abads
    end
    
    % update params
    params.mass_hip = mass_hip;
    params.mass_knee = mass_knee;
    params.mass_body = mass_body;
    params.CoM_body = CoM_body;
    params.I_hip = I_hip;
    params.I_knee = I_knee;
    params.I_body = I_body;
end

%% model initilization
robot.NB = 7;                                  % number of moving bodies (2 fictitious bodies for trunk translations)
robot.parent  = zeros(1,robot.NB);             % parent body indices
robot.Xtree   = repmat({eye(6)},robot.NB,1);   % coordinate transforms
robot.jtype   = repmat({'  '},robot.NB,1);     % joint types
robot.I       = repmat({zeros(6)},robot.NB,1); % spatial inertias
robot.gravity = [0 0 -params.g]';              % gravity acceleration vec

nb = 0;                                        % current body index

%% trunk translation x (body num 1) (massless)
nb = nb + 1;
robot.parent(nb) = nb - 1;
robot.Xtree{nb} = eye(1); 
robot.jtype{nb} = 'Px'; 
robot.I{nb} = mcI(0, zeros(3,1), zeros(3,3)); 

%% trunk translation z direction (body num 2) (massless)
nb = nb + 1;
robot.parent(nb) = nb - 1;
robot.Xtree{nb} = eye(1); 
robot.jtype{nb} = 'Pz'; 
robot.I{nb} = mcI(0, zeros(3,1), zeros(3,3)); 

%% trunck rotation about y direction (body num 3)
nb = nb + 1;
robot.parent(nb) = nb - 1;
robot.Xtree{nb} = eye(1); 
robot.jtype{nb} = 'Ry';  
robot.I{nb} = mcI( mass_body, CoM_body, I_body);

nbase = nb; % floating base index for attaching two children links (hip links)

for i = 1:NLEGS
    %% Hip link
    nb = nb + 1;
    robot.parent(nb) = nbase; % parent of the hip link is base
    robot.Xtree{nb} = plux( ry(0), [hip_x(i), 0, 0]');  % translation (half the body length)
    robot.jtype{nb} = 'Ry';
    robot.I{nb} = mcI( mass_hip, CoM_hip, I_hip); 

    %% Knee link
    nb = nb + 1;
    robot.parent(nb) = nb - 1; % parent of the knee link is hip link
    robot.Xtree{nb} = plux( ry(0),[0 0 -hiplinkLength]);    % translation (length of hip link)
    robot.jtype{nb} = 'Ry';
    robot.I{nb} = mcI( mass_knee, CoM_knee, I_knee); 
end
end
