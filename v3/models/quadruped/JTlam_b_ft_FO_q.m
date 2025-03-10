function JTlam_back_ft_FO_q = JTlam_b_ft_FO_q(q1,q2,q3,q4,q5,q6,q7,lam1,lam2)
%JTLAM_B_FT_FO_Q
%    JTLAM_BACK_FT_FO_Q = JTLAM_B_FT_FO_Q(Q1,Q2,Q3,Q4,Q5,Q6,Q7,LAM1,LAM2)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    19-Jan-2023 17:11:47

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
t30 = -lam2.*(t23-t25);
t29 = lam1.*t27;
t31 = t16+t17+t27;
t32 = t15+t20+t28;
t33 = lam1.*t31;
t34 = lam2.*t32;
t35 = t29+t30;
t36 = t33+t34;
JTlam_back_ft_FO_q = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,lam1.*(t17+t27+t2.*t21)+lam2.*(t15+t28-t5.*t21),0.0,0.0,t36,t35,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t36,0.0,0.0,t36,t35,0.0,0.0,t35,0.0,0.0,t35,t35],[7,7]);
