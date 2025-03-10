function s_ring_v_deriv = s_ring_v_deriv(in1,in2)
%S_RING_V_DERIV
%    S_RING_V_DERIV = S_RING_V_DERIV(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    31-Dec-2024 15:58:41

q2 = in1(2,:);
q3 = in1(3,:);
v1 = in2(1,:);
v2 = in2(2,:);
v3 = in2(3,:);
t2 = cos(q2);
t3 = cos(q3);
t4 = sin(q2);
t5 = sin(q3);
s_ring_v_deriv = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t4.*v1.*v2,-v1.*(t2.*t5.*v2+t3.*t4.*v3),-v1.*(t2.*t3.*v2-t4.*t5.*v3),0.0,0.0,0.0,0.0,-v1.*(t3.*t4.*v2+t2.*t5.*v3)-t3.*v2.*v3,-v1.*(t2.*t3.*v3-t4.*t5.*v2)+t5.*v2.*v3],[6,3]);
