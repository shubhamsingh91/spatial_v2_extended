
clear all; clc;

% modID derivs with a floating base
N = 17;

% Create a random model with N links
model = autoTree(N, 3, pi/3);
model.jtype{1} = 'Fb';
model.jtype{2} = 'Fb';

model = autoTree(N, 1);

model = postProcessModel(model);

% Random configuration and velocity
q   = ones(model.NQ,1);
q   = normalizeConfVec(model, q); 

qd  = ones(model.NV,1);
qdd = ones(model.NV,1);
lambda = ones(model.NV,1);
mu = ones(model.NV,1);
newConfig = @(x) configurationAddition(model,q,x);

%% Mod ID
[tau]      = ID(model, q ,qd ,qdd);                    % Inverse dynamics
out = modID( model, q, qd, qdd, lambda );              % modified inverse Dynamics

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
derivs_ground = modID_second_derivatives_ground( model, q, qd, qdd, lambda);
derivs = modID_second_derivatives( model, q, qd, qdd, lambda);
 
ID_SO_q_ground = derivs_ground.dmod_dqq;
ID_SO_v_ground = derivs_ground.dmod_dvv;
ID_SO_qv_ground = derivs_ground.dmod_dqv;
ID_SO_vq_ground = ID_SO_qv_ground.';
% ID_SO_aq_ground = derivs_ground.dmod_daq;

ID_SO_q = derivs.dmod_dqq;
ID_SO_v = derivs.dmod_dvv;
ID_SO_qv = derivs.dmod_dqv;
ID_SO_vq = ID_SO_qv.';

% complex-step
modID_cs_qq  = complexStepJacobian( @(x) outputSelect(1,@modID_derivatives_ground,...
            model,x,qd,qdd,lambda),q );
    
modID_cs_vv = complexStepJacobian( @(x) outputSelect(2,@modID_derivatives_ground,...
                model,q,x,qdd,lambda),qd );

modID_cs_qv  = complexStepJacobian( @(x) outputSelect(2,@modID_derivatives_ground,...
            model,x,qd,qdd,lambda),q );

modID_cs_vq = complexStepJacobian( @(x) outputSelect(1,@modID_derivatives_ground,...
                model,q,x,qdd,lambda),qd );

modID_cs_qa  = complexStepJacobian( @(x) outputSelect(3,@modID_derivatives_ground,...
            model,x,qd,qdd,lambda),q );
 
modID_cs_aq  = complexStepJacobian( @(x) outputSelect(1,@modID_derivatives_ground,...
            model,q,qd,x,lambda),qdd );
   
        
checkValue('modID_qq'   , ID_SO_q      , modID_cs_qq            ); % Partials of modID w.r.t. q
checkValue('modID_vv'   , ID_SO_v      , modID_cs_vv            ); % Partials of modID w.r.t. v
checkValue('modID_qv'   , ID_SO_qv      , modID_cs_qv            ); % Partials of modID w.r.t. q,v
checkValue('modID_vq'   , ID_SO_vq      , modID_cs_vq           ); % Partials of modID w.r.t. v,q
checkValue('modID_aq sanity'   , modID_cs_aq      , modID_cs_qa.'           ); % Partials of modID w.r.t. a,q

checkValue('modID_qq_ground'   , ID_SO_q_ground      , modID_cs_qq            ); % Partials of modID w.r.t. q
checkValue('modID_vv_ground'   , ID_SO_v_ground      , modID_cs_vv            ); % Partials of modID w.r.t. v
checkValue('modID_qv_ground'   , ID_SO_qv_ground      , modID_cs_qv            ); % Partials of modID w.r.t. q,v
checkValue('modID_vq_ground'   , ID_SO_vq_ground      , modID_cs_vq            ); % Partials of modID w.r.t. v,q
% checkValue('modID_aq_ground'   , ID_SO_aq_ground      , modID_cs_aq            ); % Partials of modID w.r.t. a,q

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