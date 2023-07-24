function a = Impact_dyn(robot,q,qd_pre,e)

robot =Jac_build(robot,q,qd_pre,zeros(robot.NV,1),zeros(2,1),robot.quad_foot);
[H,~] = HandC(robot,q,qd_pre);
J = robot.J;
 

K_inv = inv([H,J.';J,zeros(size(J,1))]);

a = K_inv*[H*qd_pre; e*J*qd_pre]; % [v+ -\lambda_c]


end