clear; clc
N = 1;

% Create a random model with N links
model = autoTree(N,1);

model.jtype{1} = 'S';
model = postProcessModel(model);

syms q [model.NQ, 1] real  % Create a symbolic vector q_sym

 S = sphericalZYXSubspace(q);
 syms m [3,1 ] real
 jac_expr1 = jacobian(S*m,q);

 % Convert the Jacobian into a MATLAB function and save it as a .m file
matlabFunction(jac_expr1, 'Vars', {q, m}, 'File', 'jacobian_function.m');

% Convert the Jacobian into a C++ function and save it as a .cpp file
matlabFunction(jac_expr1, 'Vars', {q, m}, 'File', 'jacobian_function_cpp', 'Optimize', true);
codegen jacobian_function -args {zeros(3, 1), zeros(3, 1)}
