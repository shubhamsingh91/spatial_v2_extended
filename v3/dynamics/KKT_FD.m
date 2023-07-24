function Xf = KKT_FD(robot,q,qd,tau,S)

[H,C] = HandC(robot,q,qd);
 robot =Jac_build(robot,q,qd,zeros(robot.NV,1),zeros(2,1),robot.quad_foot);

J = robot.J;
Jd = robot.Jd;

K_inv = inv([H,J.';J,zeros(size(J,1))]);
Xf = K_inv*[S*tau-C;-Jd*qd];                          % [qdd, -\lambda_i]

end