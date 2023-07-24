function [derivs_SO]= KKT_SO_derivatives(robot,q,qd,tau,S)
  % Tailored for the Quadruped model- Needs to be changed for different
  % model
  
 % KKT Dynamics
    robot =Jac_build(robot,q,qd,zeros(robot.NV,1),zeros(2,1),robot.quad_foot);
    a = KKT_FD(robot,q,qd,tau,S);
    qdd = a(1:robot.NV); lam=-a(robot.NV+1:end);
    robot =Jac_build(robot,q,qd,qdd,lam,robot.quad_foot);
   
    [dtau_dq,dtau_dv] = ID_derivatives(robot,q,qd,qdd);
    
    derivs = ID_SO_derivatives( robot, q, qd, qdd);
    
    d2tau_dq = derivs.d2tau_dq;
    d2tau_dv = derivs.d2tau_dv;
    d2tau_dqv = derivs.d2tau_dqv;
    M_FO = derivs.dM_dq;
    d2tau_dvq = rotR(d2tau_dqv);
    
    [H,C] = HandC(robot,q,qd);
    J = robot.J;
    Jd = robot.Jd;
    K_inv = inv([H,J.';J,zeros(size(J,1))]);

    dKKT_dq = -K_inv*[dtau_dq-robot.JTlam_FO;...
                         robot.Jqdd_FO_q+robot.Jdqd_FO_q]; 

    dKKT_dv = -K_inv*[dtau_dv;robot.Jdqd_FO_v];

    dKKT_dtau = K_inv*[S;zeros(2,4)];
    
     FD_FO_q =        dKKT_dq(1:robot.NV,:);
     FD_FO_lam_q =    -dKKT_dq(robot.NV+1:end,:);
     FD_FO_v =        dKKT_dv(1:robot.NV,:);
     FD_FO_lam_v =    -dKKT_dv(robot.NV+1:end,:);
     FD_FO_tau =      dKKT_dtau(1:robot.NV,:);
     FD_FO_lam_tau =  -dKKT_dtau(robot.NV+1:end,:);
 
   term1 = [d2tau_dq-robot.JTlam_SOq;robot.Jqdd_SO_q+robot.Jdqd_SO_q];
   
   term2 = [Tm(M_FO,FD_FO_q)+rotR(Tm(M_FO,FD_FO_q))-...
                Tm(robot.JT_FO_q,FD_FO_lam_q)-rotR(Tm(robot.JT_FO_q,FD_FO_lam_q));...
        Tm(robot.J_FO_q,FD_FO_q)+rotR(Tm(robot.J_FO_q,FD_FO_q))];
    
  derivs_SO.da_dqq = -mT(K_inv,term1+term2);
  derivs_SO.da_dvv = -mT(K_inv,[d2tau_dv;robot.Jdqd_SO_v]);
  derivs_SO.da_dvq = -mT(K_inv,[d2tau_dvq+Tm(M_FO,FD_FO_v)- ...
                    Tm(robot.JT_FO_q,FD_FO_lam_v);robot.Jdqd_SO_vq + Tm(robot.J_FO_q,FD_FO_v)]);
  derivs_SO.da_dqv =rotR(derivs_SO.da_dvq);
  derivs_SO.da_dtauq = -mT(K_inv,[Tm(M_FO,FD_FO_tau) - Tm(robot.JT_FO_q,FD_FO_lam_tau)...
                                    ; Tm(robot.J_FO_q,FD_FO_tau)]);
                                
  derivs_SO.da_dqtau = rotR(derivs_SO.da_dtauq);

end
