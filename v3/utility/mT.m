function tens = mT(m,T)
% Description- 
% Product of a matrix and a tensor
% Inputs 
% m (N1XN2), T(N2XN3XN4)
% Outputs
% tens (N1XN3XN4)

for ii=1:size(T,3)
    tens(:,:,ii) = m*T(:,:,ii);
end