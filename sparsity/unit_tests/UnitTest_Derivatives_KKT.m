
% Unit tests for derivatives functions of KKT & Impact Dynamics

clear

% Create a random model with N links
[robot, robot_params] = build_quadruped_model('2D');
robot = postProcessModel(robot);

checkDerivatives(robot,robot_params,'Floating Base No Rotors');

function checkDerivatives(robot,robot_params, desc)

    % Random configuration and velocity
    u_size=4;
    q   = rand(robot.NQ,1);
    qd  = rand(robot.NV,1);
    tau = rand(u_size,1);
    
    
    %  selection matrix
    S = [zeros(robot.NQ-u_size,...
    u_size);
    eye(u_size)];
 
    robot.quad_foot=1; % front=1, back=2
  
    % Check KKT dynamics FO derivatives
   [da_dq,da_dv,da_dtau]=KKT_FO_derivatives(robot,q,qd,tau,S);
   
    da_dq_cs  = complexStepJacobian( @(y) KKT_FD(robot,y,qd,tau,S),q);
    da_dv_cs  = complexStepJacobian( @(y) KKT_FD(robot,q,y,tau,S),qd);
    da_dtau_cs  = complexStepJacobian( @(y) KKT_FD(robot,q ,qd,y,S),tau);
      
    % Check KKT dynamics SO derivatives
   [derivs_KKT]= KKT_SO_derivatives(robot,q,qd,tau,S);
 
   [da_dqq_cs, da_dvq_cs, da_dtauq_cs] = complexStepHessian_FD(@(x) KKT_FO_derivatives(robot, x ,qd ,tau,S) , q);
   [da_dqv_cs, da_dvv_cs, da_dtauv_cs] = complexStepHessian_FD(@(x) KKT_FO_derivatives(robot, q ,x ,tau,S) , qd);
   [da_dqtau_cs, da_dvtau_cs, da_dtautau_cs] = complexStepHessian_FD(@(x) KKT_FO_derivatives(robot, q ,qd ,x,S) , tau);

    % Check Impact Dynamics FO Derivatives
    [d_imp_dq,d_imp_dv]=Impact_FO_derivatives(robot,q,qd);
   
  d_imp_dq_cs  = complexStepJacobian( @(y) Impact_dyn(robot,y,qd,0),q);
  d_imp_dv_cs  = complexStepJacobian( @(y) Impact_dyn(robot,q,y,0),qd);
              
   % Check Impact Dynamics SO Derivatives
   derivs_impact = Impact_SO_derivatives(robot,q,qd);
   
  [d_imp_dqq_cs, d_imp_dvq_cs] = complexStepHessian(@(x) Impact_FO_derivatives(robot, x ,qd) , q);    
  [d_imp_dqv_cs, d_imp_dvv_cs] = complexStepHessian(@(x) Impact_FO_derivatives(robot, q ,x) , qd);

 checkValue('KKT_q'   , da_dq      , da_dq_cs   );                                            % Partials of FD w.r.t. q
 checkValue('KKT_v'  ,  da_dv     , da_dv_cs  );                                              % Partials of FD w.r.t. v
 checkValue('KKT_tau' , da_dtau    , da_dtau_cs );                                            % Partials of FD w.r.t. tau

 checkValue('KKT_qq'       , derivs_KKT.da_dqq         , da_dqq_cs                );          % SO Partials of FD w.r.t. q,q
 checkValue('KKT_vv'       , derivs_KKT.da_dvv         , da_dvv_cs                );          % SO Partials of FD w.r.t. v,v
 checkValue('KKT_qv'       , derivs_KKT.da_dqv         , da_dqv_cs                );          % SO Partials of FD w.r.t. q,v
 checkValue('KKT_vq'       , derivs_KKT.da_dvq         , da_dvq_cs                );          % SO Partials of FD w.r.t. v,q
 checkValue('KKT_qtau'       , derivs_KKT.da_dtauq         , da_dtauq_cs                );    % SO Partials of FD w.r.t. tau,q
 checkValue('KKT_tauq'       , derivs_KKT.da_dqtau         , da_dqtau_cs                );    % SO Partials of FD w.r.t. q,tau

 checkValue('Impact_q'   , d_imp_dq      , d_imp_dq_cs   );                                   % Partials of Impact w.r.t. q
 checkValue('Impact_v'  ,  d_imp_dv     , d_imp_dv_cs  );                                     % Partials of Impact w.r.t. v

 checkValue('Impact_qq'       , derivs_impact.da_dqq         , d_imp_dqq_cs                ); % SO Partials of Impact w.r.t. q,q
 checkValue('Impact_vv'       , derivs_impact.da_dvv         , d_imp_dvv_cs                ); % SO Partials of Impact w.r.t. v,v
 checkValue('Impact_qv'       , derivs_impact.da_dqv         , d_imp_dqv_cs                ); % SO Partials of Impact w.r.t. q,v
 checkValue('Impact_vq'       , derivs_impact.da_dvq         , d_imp_dvq_cs                ); % SO Partials of Impact w.r.t. v,q

  fprintf('\n');
end



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