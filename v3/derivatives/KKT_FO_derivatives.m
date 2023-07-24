function [dKKT_dq,dKKT_dv,dKKT_dtau]= KKT_FO_derivatives(robot,q,qd,tau,S)
  % Tailored for the Quadruped model- Needs to be changed for different
  % model
  
 % KKT Dynamics
    robot =Jac_build(robot,q,qd,zeros(robot.NV,1),zeros(2,1),robot.quad_foot);
    a = KKT_FD(robot,q,qd,tau,S);
    qdd = a(1:robot.NV); lam=-a(robot.NV+1:end);
    robot =Jac_build(robot,q,qd,qdd,lam,robot.quad_foot);
    [dtau_dq,dtau_dv] = ID_derivatives(robot,q,qd,qdd);

    [H,C] = HandC(robot,q,qd);
    J = robot.J;
    Jd = robot.Jd;
    K_inv = inv([H,J.';J,zeros(size(J,1))]);

    dKKT_dq = -K_inv*[dtau_dq-robot.JTlam_FO;...
                         robot.Jqdd_FO_q+robot.Jdqd_FO_q]; 

    dKKT_dv = -K_inv*[dtau_dv;robot.Jdqd_FO_v];

    dKKT_dtau = K_inv*[S;zeros(2,4)];

end
