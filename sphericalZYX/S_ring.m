function s_ring = s_ring(in1,in2)
%S_RING
%    S_RING = S_RING(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    30-Dec-2024 18:45:44

q2 = in1(2,:);
q3 = in1(3,:);
v2 = in2(2,:);
v3 = in2(3,:);
t2 = cos(q2);
t3 = cos(q3);
t4 = sin(q2);
t5 = sin(q3);
s_ring = reshape([0.0,0.0,0.0,-t2.*v2,t2.*t3.*v3-t4.*t5.*v2,-t3.*t4.*v2-t2.*t5.*v3,0.0,0.0,0.0,0.0,-t5.*v3,-t3.*v3,0.0,0.0,0.0,0.0,0.0,0.0],[6,3]);
