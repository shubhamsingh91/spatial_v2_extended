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
syms f [6,1 ] real
syms m [3,1 ] real
dSm_sq = jacobian(S*m,q); 

%% S_ring

for i = 1:3
    temp1 = jacobian(S(:,i),q); % this tensor is d(S_i)/dq , each page is d(S_i)/dq
    s_ring(:,i) = temp1*v;% this is d(S_i)/dq  * v
end

% Pat's version
s_ring_v2 = sym(zeros(6,3));
for i = 1:3
    ds_dq(:,:,i) = diff(S,q(i)); % this tensor is d(S_i)/dq , each page is d(S_i)/dq
    s_ring_v2  =s_ring_v2+ ds_dq(:,:,i)*v(i);% this is d(S_i)/dq  * v
end

s_ring_v3 = jacobian(jacobian(S*m,q)*v,m);

disp(s_ring - s_ring_v2)
disp(s_ring - s_ring_v3)

%sanity check
disp("sanity check = ")
disp(jacobian(s_ring*v,v) - (s_ring + jacobian(S*v,q)))

matlabFunction(s_ring, 'Vars', {q, v}, 'File', 's_ring', 'Outputs', {'s_ring'});
%

%% d/dq Sm, d/dq S'T * f

dSm_dq = jacobian(S*m,q);
matlabFunction(dSm_dq, 'Vars', {q, m}, 'File', 'dSm_dq', 'Outputs', {'dSm_dq'});

dSTf_dq = jacobian(S.'*f,q);
matlabFunction(dSTf_dq, 'Vars', {q, f}, 'File', 'dSTf_dq', 'Outputs', {'dSTf_dq'});

%% partial(S_ring * v)/partial q

s_ring_v_deriv = jacobian(s_ring*v,q);
matlabFunction(s_ring_v_deriv, 'Vars', {q, v}, 'File', 's_ring_v_deriv', 'Outputs', {'s_ring_v_deriv'});


%% testing numerical values
qnum = [1,1,1].'; mnum = [1,2,3].';
vnum = [1,1,1].'; fnum = [1,2,3,4,5,6].';

s_ring_num  = double(subs(s_ring,[q,v],[qnum,vnum]))

s_ring_v_deriv_num = double(subs(s_ring_v_deriv,[q,v],[qnum,vnum]))

dsm_dq_num = double(subs(dSm_dq,[q,m],[qnum,mnum]))
