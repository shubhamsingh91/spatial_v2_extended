clear; clc
N = 1;

% Create a random model with N links
model = autoTree(N,1);

model.jtype{1} = 'S';
model = postProcessModel(model);

syms q [3, 1] real  % Create a symbolic vector q_sym
S = sphericalZYXSubspace(q);

% S * m
syms m [3,1 ] real
jac_expr1 = jacobian(S*m,q);
matlabFunction(jac_expr1, 'Vars', {q, m}, 'File', 'S_ring', 'Outputs', {'jac_expr1'}, 'Vars', {q, m});

% S.' * m
syms n [6,1 ] real
jac_expr2 = jacobian(S.'*n,q);
matlabFunction(jac_expr2, 'Vars', {q, n}, 'File', 'STn.m');


