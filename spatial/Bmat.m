function B = Bmat(I,v)

% Generalized Coriolis Matrix
% Inputs - 
% I (6x6 Inertia matrix), v (6x1 spatial motion vector)
% Outputs-
% B (6x6 Coriolis matrix)

B = 0.5*(crf(v)*I-I*(crm(v))+cmf_bar(I*v)); 


end