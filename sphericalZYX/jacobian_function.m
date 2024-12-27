function jac_expr1 = jacobian_function(in1, in2)
%JACOBIAN_FUNCTION
%    JAC_EXPR1 = JACOBIAN_FUNCTION(IN1,IN2)

%    This function computes the Jacobian based on inputs in1 and in2.

% Inputs:
%   in1: 2x1 vector of generalized coordinates (q1, q2)
%   in2: 3x1 vector of some parameters (m1, m2, m3)

% Output:
%   jac_expr1: 6x4 matrix representing the Jacobian

m1 = in2(1,:);  % m1 value from in2
m2 = in2(2,:);  % m2 value from in2
m3 = in2(3,:);  % m3 value from in2
q1 = in1(1,:);  % q1 value from in1
q2 = in1(2,:);  % q2 value from in1

% Pre-compute the sines and cosines
t2 = cos(q1);  % cos(q1)
t3 = cos(q2);  % cos(q2)
t4 = sin(q1);  % sin(q1)
t5 = sin(q2);  % sin(q2)

% Calculate the Jacobian matrix
jac_expr1 = reshape([0.0, 0.0, 0.0, -m1.*t2 - m3.*t3.*t4 - m2.*t4.*t5, ...
                     0.0, 0.0, 0.0, 0.0, 0.0, ...
                     m2.*t2.*t3 - m3.*t2.*t5, -m3.*t3 - m2.*t5, ...
                     0.0, 0.0, 0.0, 0.0, 0.0, ...
                     0.0, 0.0, 0.0, 0.0, 0.0, ...
                     0.0, 0.0, 0.0, 0.0], [6, 4]);
end
