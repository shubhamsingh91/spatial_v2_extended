function Jdqd_front_ft_SO_v = Jdqd_fr_ft_SO_v(q1,q2,q3,q4,q5,q6,q7,qd1,qd2,qd3,qd4,qd5,qd6,qd7)
%JDQD_FR_FT_SO_V
%    JDQD_FRONT_FT_SO_V = JDQD_FR_FT_SO_V(Q1,Q2,Q3,Q4,Q5,Q6,Q7,QD1,QD2,QD3,QD4,QD5,QD6,QD7)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    19-Jan-2023 17:11:43

t2 = cos(q3);
t3 = cos(q4);
t4 = cos(q5);
t5 = sin(q3);
t6 = sin(q4);
t7 = sin(q5);
t8 = t3.*t4;
t9 = t3.*t7;
t10 = t4.*t6;
t11 = t6.*t7;
t13 = t6.*(2.09e+2./1.0e+3);
t15 = t2.*t3.*(2.09e+2./5.0e+2);
t16 = t2.*t6.*(2.09e+2./5.0e+2);
t17 = t3.*t5.*(2.09e+2./5.0e+2);
t18 = t5.*t6.*(2.09e+2./5.0e+2);
t12 = -t11;
t14 = t9+t10;
t19 = -t18;
t21 = t13-1.9e+1./1.0e+2;
t20 = t8+t12;
t22 = t2.*t14.*(3.9e+1./1.0e+2);
t23 = t5.*t14.*(3.9e+1./1.0e+2);
t24 = -t23;
t25 = t2.*t20.*(3.9e+1./1.0e+2);
t26 = t5.*t20.*(3.9e+1./1.0e+2);
t27 = t22+t26;
t28 = t24+t25;
t29 = t16+t17+t27;
t30 = t15+t19+t28;
Jdqd_front_ft_SO_v = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t17+t27+t2.*t21.*2.0,t15+t28-t5.*t21.*2.0,t29,t30,t27,t28,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t29,t30,t29,t30,t27,t28,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t27,t28,t27,t28,t27,t28,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[2,7,7]);
