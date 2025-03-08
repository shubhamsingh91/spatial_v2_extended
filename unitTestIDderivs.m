
clear all; clc;

% modID derivs with a floating base
N = 17;

% Create a random model with N links
model = autoTree(N, 3, pi/3);
model.jtype{1} = 'Fb';
model.jtype{2} = 'Fb';

model = postProcessModel(model);

% Random configuration and velocity
q   = rand(model.NQ,1);
q   = normalizeConfVec(model, q); 

qd  = rand(model.NV,1);
qdd = rand(model.NV,1);
lambda = rand(model.NV,1);
newConfig = @(x) configurationAddition(model,q,x);

%% Mod ID Derivs

[dmodID_dq, dmodID_dqd] = modID_derivatives( model, q, qd, qdd, lambda );

dmodID_dq_cs  = complexStepJacobian(@(x) modID(model, newConfig(x) , ...
                qd ,qdd,lambda), zeros(model.NV,1) );
dmodID_dqd_cs = complexStepJacobian(@(x) modID(model, q ,x  ,qdd,lambda), qd);

checkValue('modID_q'   , dmodID_dq      , dmodID_dq_cs            ); % Partials of modID w.r.t. q
checkValue('modID_qd'  , dmodID_dqd     , dmodID_dqd_cs           ); % Partials of modID w.r.t. qd

[tau]      = ID(model, q ,qd ,qdd);                    % Inverse dynamics

%% Mod FD derivs
[dmodFD_dq, dmodFD_dqd, dmodFD_dtau] = modFD_derivatives( model, q, qd, tau, lambda );

dmodFD_dq_cs  = complexStepJacobian(@(x) modFD(model,  newConfig(x) ,qd ,tau,lambda), zeros(model.NV,1) );
dmodFD_dqd_cs = complexStepJacobian(@(x) modFD(model, q ,x  ,tau,lambda), qd);
dmodFD_dtau_cs = complexStepJacobian(@(x) modFD(model,q ,qd ,x ,lambda), tau);


checkValue('modFD_q'   , dmodFD_dq      , dmodFD_dq_cs            ); % Partials of modFD w.r.t. q
checkValue('modFD_qd'  , dmodFD_dqd     , dmodFD_dqd_cs           ); % Partials of modFD w.r.t. qd
checkValue('modFD_tau' , dmodFD_dtau    , dmodFD_dtau_cs          ); % Partials of modFD w.r.t. qd

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