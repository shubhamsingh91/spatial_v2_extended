function [derivs_SO]= Impact_SO_derivatives(robot,q,qd_pre)
  % Tailored for the Quadruped model- Needs to be changed for different
  % model
  
    qdd = Impact_dyn(robot,q,qd_pre,0);      % [v+ -\lambda_c] 
    qd_post = qdd(1:robot.NV);               % v after the impact
    lam = -qdd(8:9);                        % [\lambda_i]
    dv = qd_post;
    
     [d_imp_dq_full,d_imp_dv_full]=Impact_FO_derivatives(robot,q,qd_pre);
     
     FD_FO_q = d_imp_dq_full(1:robot.NV,:);
     FD_FO_lam_q = -d_imp_dq_full(robot.NV+1:end,:);
     FD_FO_v =        d_imp_dv_full(1:robot.NV,:);
     FD_FO_lam_v =    -d_imp_dv_full(robot.NV+1:end,:);
     
    [H,~] = HandC(robot,q,qd_pre);
    
    robot =Jac_build(robot,q,qd_pre,dv,lam,robot.quad_foot);
    robot_no_grav = robot;
    robot_no_grav.gravity = zeros(3,1);
    J = robot.J;
    
    derivs = ID_SO_derivatives( robot_no_grav, q, zeros(robot.NV,1), qd_post-qd_pre);
    term2_1 = Tm(derivs.dM_dq,FD_FO_q);
    term1 = [derivs.d2tau_dq-robot.JTlam_SOq;robot.Jqdd_SO_q];
  
    term2 = [term2_1+rotR(term2_1)-...
                    Tm(robot.JT_FO_q,FD_FO_lam_q)-rotR(Tm(robot.JT_FO_q,FD_FO_lam_q));...
            Tm(robot.J_FO_q,FD_FO_q)+rotR(Tm(robot.J_FO_q,FD_FO_q))];
        
   K_inv = inv([H,J.';J,zeros(size(J,1))]);
     
  derivs_SO.da_dqq = -mT(K_inv,term1+term2);
  derivs_SO.da_dvv = zeros(robot.NV+2,robot.NV,robot.NV);
  derivs_SO.da_dvq = mT(K_inv,[Tm(derivs.dM_dq,eye(robot.NV))-Tm(derivs.dM_dq,FD_FO_v)+ ...
                Tm(robot.JT_FO_q,FD_FO_lam_v);0*Tm(robot.J_FO_q,eye(robot.NV))-Tm(robot.J_FO_q,FD_FO_v)]);
  derivs_SO.da_dqv = rotR(derivs_SO.da_dvq);


end
