function tens = Tm(T,m)
% Description- 
% Product of a tensor and a matrix
% Inputs 
% T(N1XN2XN3), m(N2XN4)
% Outputs
% tens (N1XN4XN3)

  for ii=1:size(T,3)  
   tens(:,:,ii) = T(:,:,ii)*m;
  end
