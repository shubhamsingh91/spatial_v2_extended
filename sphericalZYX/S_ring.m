function jac_expr1 = S_ring(q,m)
%S_RING
%    JAC_EXPR1 = S_RING(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Dec-2024 11:33:20

m1 = m(1,:);
m2 = m(2,:);
q2 = q(2,:);
q3 = q(3,:);
t2 = cos(q2);
t3 = cos(q3);
t4 = sin(q2);
t5 = sin(q3);
jac_expr1 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-m1.*t2,-m1.*t4.*t5,-m1.*t3.*t4,0.0,0.0,0.0,0.0,-m2.*t5+m1.*t2.*t3,-m2.*t3-m1.*t2.*t5],[6,3]);