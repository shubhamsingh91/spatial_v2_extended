function tens_out = rotT(tens_in)

% Description- 
% rotT operator on a tensor
%(tens_in) \Tten
% Inputs 
% tens_in (N1XN2XN3)
% Outputs
% tens_out(N2XN1XN3 tensor)

for ii=1:size(tens_in,3)
    tens_out(:,:,ii) = tens_in(:,:,ii).';
end