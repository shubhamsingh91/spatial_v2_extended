function [da_dq,da_dv]= Impact_FO_derivatives(robot,q,v_pre)

 % Impact Dynamics- Assuming coeff. of restitution = 0

    [H,C] = HandC(robot,q,v_pre);
    
    qdd = Impact_dyn(robot,q,v_pre,0);      % [v+ -\lambda_c] 
    v_post = qdd(1:robot.NB);               % v after the impact
    lam = -qdd(8:9);                        % [\lambda_i]
    dv = v_post;
    
    robot =Jac_build(robot,q,dv,dv,lam,robot.quad_foot);
    
    J = robot.J;
    K_inv = inv([H,J.';J,zeros(size(J,1))]);
    
    robot_no_grav = robot;
    robot_no_grav.gravity = zeros(3,1);
    
    [dtau_dq, dtau_dqd] = ID_derivatives(robot_no_grav, q, zeros(robot.NV,1), v_post - v_pre ); % IDFOZA 
   
    da_dq = -K_inv*[dtau_dq-robot.JTlam_FO;robot.Jqdd_FO_q]; 
    da_dv = K_inv*[H; 0*J];


end
