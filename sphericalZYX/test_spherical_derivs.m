clear; clc
N = 1;

% Create a random model with N links
model = autoTree(N,1);

model.jtype{1} = 'S';
model = postProcessModel(model);

syms q [3, 1] real  % Create a symbolic vector q_sym
S = sphericalZYXSubspace(q);

% S_ring = partial (S)/partial (q)  * v
syms v [3,1 ] real
s_ring_jac = jacobian(S*v,q); % not correct

for i = 1:3
    temp1 = jacobian(S(:,i),q); % this tensor is d(S_i)/dq , each page is d(S_i)/dq
    s_ring(:,i) = temp1*v;% this is d(S_i)/dq  * v
end

% matlabFunction(s_ring, 'Vars', {q, v}, 'File', 's_ring', 'Outputs', {'s_ring'});
% 
% partial (S_ring * v)/partial q
% s_ring_times_v_partial = jacobian(s_ring_v2*v,q);
% matlabFunction(s_ring_times_v_partial, 'Vars', {q, v}, 'File', 's_ring_v_deriv', 'Outputs', {'s_ring_v_partial'});


