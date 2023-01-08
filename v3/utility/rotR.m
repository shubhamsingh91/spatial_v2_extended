function tens_out = rotR(tens_in)

% Description- 
% rotR operator on a tensor
%(tens_in) \Rten
% Inputs 
% tens_in (N1XN2XN3)
% Outputs
% tens_out(N1XN3XN2 tensor)

m = size(tens_in,3); % pages
n = size(tens_in,2); % cols

 for jj=1:m
   for ii=1:n
        temp1 = tens_in(:,ii,jj);
        tens_out(:,jj,ii) = temp1;
    end    
end