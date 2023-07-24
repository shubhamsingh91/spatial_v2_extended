function Jqdd_back_ft_FO_q = Jqdd_b_ft_FO_q(q1,q2,q3,q4,q5,q6,q7,qdd1,qdd2,qdd3,qdd4,qdd5,qdd6,qdd7)
%JQDD_B_FT_FO_Q
%    JQDD_BACK_FT_FO_Q = JQDD_B_FT_FO_Q(Q1,Q2,Q3,Q4,Q5,Q6,Q7,QDD1,QDD2,QDD3,QDD4,QDD5,QDD6,QDD7)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    19-Jan-2023 17:11:49

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
t30 = -qdd7.*(t23-t25);
t29 = qdd7.*t27;
t31 = t16+t17+t27;
t32 = t15+t20+t28;
t33 = qdd6.*t31;
t34 = qdd6.*t32;
Jqdd_back_ft_FO_q = reshape([0.0,0.0,0.0,0.0,t29+t33+qdd3.*(t17+t27+t2.*t21),t30+t34+qdd3.*(t15+t28-t5.*t21),0.0,0.0,0.0,0.0,t29+t33+qdd3.*t31,t30+t34+qdd3.*t32,t29+qdd3.*t27+qdd6.*t27,t30-qdd3.*(t23-t25)-qdd6.*(t23-t25)],[2,7]);
