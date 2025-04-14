
clear all; clc;

% modID derivs with a floating base
N = 17;

% Create a random model with N links
model = autoTree(N, 3, pi/3);
model.jtype{1} = 'Fb';
model.jtype{2} = 'Fb';

model = postProcessModel(model);

% Random configuration and velocity
q   = ones(model.NQ,1);
q   = normalizeConfVec(model, q); 

qd  = ones(model.NV,1);
qdd = ones(model.NV,1);
lambda = ones(model.NV,1);
mu = rand(model.NV,1);
newConfig = @(x) configurationAddition(model,q,x);

%% Mod ID
[tau]      = ID(model, q ,qd ,qdd);                    % Inverse dynamics
out = modID( model, q, qd, qdd, lambda );               % modified inverse Dynamics

checkValue('modID'   , out      , lambda.'*tau            ); % modID

% mod FD
out = modFD( model, q, qd, tau, mu );
qdd = FDab(model,q,qd,tau);

checkValue('modFD'   , out      , mu.'*qdd            ); % modFD

%% Mod ID Derivs- FO
[ID_q, ID_v] = ID_derivatives( model, q, qd, qdd );
ID_a = CRBA(model,q);

[dmodID_dq, dmodID_dqd] = modID_derivatives( model, q, qd, qdd, lambda );
[dmodID_dq_ground, dmodID_dqd_ground, dmodID_dqdd_ground] = modID_derivatives_ground( model, q, qd, qdd, lambda );

dmodID_dq_cs  = complexStepJacobian(@(x) modID(model, newConfig(x) , ...
                qd ,qdd,lambda), zeros(model.NV,1) );
dmodID_dqd_cs = complexStepJacobian(@(x) modID(model, q ,x  ,qdd,lambda), qd);
dmodID_dqdd_cs = complexStepJacobian(@(x) modID(model, q ,qd  ,x,lambda), qdd);

checkValue('modID_q'   , dmodID_dq      , dmodID_dq_cs            ); % Partials of modID w.r.t. q
checkValue('modID_qd'  , dmodID_dqd     , dmodID_dqd_cs           ); % Partials of modID w.r.t. qd
checkValue('modID_q_ground'   , dmodID_dq_ground      , dmodID_dq_cs            ); % Partials of modID w.r.t. q
checkValue('modID_qd_ground'   , dmodID_dqd_ground      , dmodID_dqd_cs            ); % Partials of modID w.r.t. qd
checkValue('modID_qdd_ground'   , dmodID_dqdd_ground      , dmodID_dqdd_cs            ); % Partials of modID w.r.t. qdd

checkValue('modID_q_ground'   , dmodID_dq_ground      , lambda.'*ID_q            ); % Partials of modID w.r.t. q
checkValue('modID_qd_ground'   , dmodID_dqd_ground      , lambda.'*ID_v            ); % Partials of modID w.r.t. qd
checkValue('modID_qdd_ground'   , dmodID_dqdd_ground      , lambda.'*ID_a           ); % Partials of modID w.r.t. qdd

%% Mod FD derivs- FO
[FD_q, FD_v, FD_tau] = FD_derivatives( model, q, qd, qdd );

[dmodFD_dq, dmodFD_dqd, dmodFD_dtau] = modFD_derivatives_ground( model, q, qd, tau, lambda );

dmodFD_dq_cs  = complexStepJacobian(@(x) modFD(model,  newConfig(x) ,qd ,tau,lambda), zeros(model.NV,1) );
dmodFD_dqd_cs = complexStepJacobian(@(x) modFD(model, q ,x  ,tau,lambda), qd);
dmodFD_dtau_cs = complexStepJacobian(@(x) modFD(model,q ,qd ,x ,lambda), tau);


checkValue('modFD_q'   , dmodFD_dq      , dmodFD_dq_cs            ); % Partials of modFD w.r.t. q
checkValue('modFD_qd'  , dmodFD_dqd     , dmodFD_dqd_cs           ); % Partials of modFD w.r.t. qd
checkValue('modFD_tau' , dmodFD_dtau    , dmodFD_dtau_cs          ); % Partials of modFD w.r.t. qd

%% Mod ID Derivs - SO
derivs = modID_second_derivatives_ground( model, q, qd, qdd, lambda);

ID_SO_q = derivs.dmod_dqq;
ID_SO_v = derivs.dmod_dvv;
ID_SO_qv = derivs.dmod_dqv;


%% MoD FD Derivs - SO

function checkValue(name, v1, v2, tolerance)
    if nargin == 3
        tolerance = sqrt(eps);
    end
    value = norm(v1(:)-v2(:));
    fprintf('%10s \t %e\n',name,value);
    if value > tolerance
        error('%s is out of tolerance',name);
    end
end