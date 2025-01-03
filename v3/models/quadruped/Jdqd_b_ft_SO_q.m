function Jdqd_back_ft_SO_q = Jdqd_b_ft_SO_q(q1,q2,q3,q4,q5,q6,q7,qd1,qd2,qd3,qd4,qd5,qd6,qd7)
%JDQD_B_FT_SO_Q
%    JDQD_BACK_FT_SO_Q = JDQD_B_FT_SO_Q(Q1,Q2,Q3,Q4,Q5,Q6,Q7,QD1,QD2,QD3,QD4,QD5,QD6,QD7)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    19-Jan-2023 17:11:45

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
t49 = -qd3.*(qd3.*(t23-t25)+qd6.*(t23-t25)+qd7.*(t23-t25));
t50 = -qd6.*(qd3.*(t23-t25)+qd6.*(t23-t25)+qd7.*(t23-t25));
t51 = -qd7.*(qd3.*(t23-t25)+qd6.*(t23-t25)+qd7.*(t23-t25));
t52 = qd3.*(qd3.*(t23-t25)+qd6.*(t23-t25)+qd7.*(t23-t25));
t53 = qd6.*(qd3.*(t23-t25)+qd6.*(t23-t25)+qd7.*(t23-t25));
t54 = qd7.*(qd3.*(t23-t25)+qd6.*(t23-t25)+qd7.*(t23-t25));
t29 = qd3.*t27;
t30 = qd6.*t27;
t31 = qd7.*t27;
t35 = t16+t17+t27;
t36 = t15+t20+t28;
t45 = t32+t33+t34;
t66 = t52+t53+t54;
t37 = qd3.*t35;
t38 = qd6.*t35;
t39 = qd3.*t36;
t40 = qd6.*t36;
t41 = t29+t30+t31;
t42 = qd3.*t41;
t43 = qd6.*t41;
t44 = qd7.*t41;
t55 = t31+t37+t38;
t60 = t34+t39+t40;
t46 = -t42;
t47 = -t43;
t48 = -t44;
t56 = qd3.*t55;
t57 = qd6.*t55;
t61 = qd3.*t60;
t62 = qd6.*t60;
t58 = -t56;
t59 = -t57;
t63 = -t61;
t64 = -t62;
t65 = t46+t47+t48;
t67 = t48+t58+t59;
t68 = t54+t63+t64;
Jdqd_back_ft_SO_q = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t48+t59-qd3.*(t31+t38+qd3.*(t17+t27+t2.*t21)),t54+t64-qd3.*(t34+t40+qd3.*(t15+t28-t5.*t21)),0.0,0.0,0.0,0.0,t67,t68,t65,t66,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t67,t68,0.0,0.0,0.0,0.0,t67,t68,t65,t66,0.0,0.0,0.0,0.0,t65,t66,0.0,0.0,0.0,0.0,t65,t66,t65,t66],[2,7,7]);
