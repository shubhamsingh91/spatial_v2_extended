function B = Bten(I,V)
% Generalized Coriolis Tensor
% Inputs - 
% I (6x6 Inertia matrix), v (6xN spatial motion matrix)
% Outputs-
% B (6x6xN Coriolis Tensor)

B = 0.5*(Tm(cmfM(V),I)-mT(I,crmM(V))+cmf_barM(I*V));

end