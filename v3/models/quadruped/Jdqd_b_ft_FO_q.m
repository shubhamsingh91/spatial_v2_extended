function Jdqd_back_ft_FO_q = Jdqd_b_ft_FO_q(q1,q2,q3,q4,q5,q6,q7,qd1,qd2,qd3,qd4,qd5,qd6,qd7)
%JDQD_B_FT_FO_Q
%    JDQD_BACK_FT_FO_Q = JDQD_B_FT_FO_Q(Q1,Q2,Q3,Q4,Q5,Q6,Q7,QD1,QD2,QD3,QD4,QD5,QD6,QD7)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    19-Jan-2023 17:11:44

t2 = cos(q3);
t3 = cos(q6);
t4 = cos(q7);
t5 = sin(q3);
t6 = sin(q6);
t7 = sin(q7);
t8 = t3.*t4;
t9 = t3.*t7;
t10 = t4.*t6;
t11 = t6.*t7;
t13 = t6.*(2.09e+2./1.0e+3);
t15 = t2.*t3.*(2.09e+2./1.0e+3);
t17 = t3.*t5.*(2.09e+2./1.0e+3);
t20 = t5.*t6.*(-2.09e+2./1.0e+3);
t12 = -t11;
t14 = t9+t10;
t16 = t2.*t13;
t18 = t5.*t13;
t21 = t13+1.9e+1./1.0e+2;
t19 = t8+t12;
t22 = t2.*t14.*(3.9e+1./2.0e+2);
t23 = t5.*t14.*(3.9e+1./2.0e+2);
t24 = -t23;
t25 = t2.*t19.*(3.9e+1./2.0e+2);
t26 = t5.*t19.*(3.9e+1./2.0e+2);
t27 = t22+t26;
t28 = t24+t25;
t32 = -qd3.*(t23-t25);
t33 = -qd6.*(t23-t25);
t34 = -qd7.*(t23-t25);
t45 = -qd7.*(qd3.*(t23-t25)+qd6.*(t23-t25)+qd7.*(t23-t25));
t29 = qd3.*t27;
t30 = qd6.*t27;
t31 = qd7.*t27;
t35 = t16+t17+t27;
t36 = t15+t20+t28;
t43 = t32+t33+t34;
t37 = qd3.*t35;
t38 = qd6.*t35;
t39 = qd3.*t36;
t40 = qd6.*t36;
t41 = t29+t30+t31;
t42 = qd7.*t41;
t46 = t31+t37+t38;
t49 = t34+t39+t40;
t44 = -t42;
t47 = qd6.*t46;
t50 = qd6.*t49;
t48 = -t47;
Jdqd_back_ft_FO_q = reshape([0.0,0.0,0.0,0.0,t45+t50+qd3.*(t34+t40+qd3.*(t15+t28-t5.*t21)),t44+t48-qd3.*(t31+t38+qd3.*(t17+t27+t2.*t21)),0.0,0.0,0.0,0.0,t45+t50+qd3.*t49,t44+t48-qd3.*t46,t45-qd3.*(qd3.*(t23-t25)+qd6.*(t23-t25)+qd7.*(t23-t25))-qd6.*(qd3.*(t23-t25)+qd6.*(t23-t25)+qd7.*(t23-t25)),t44-qd3.*t41-qd6.*t41],[2,7]);
