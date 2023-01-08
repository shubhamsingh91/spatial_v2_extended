function mat_out = hess_to_mat2(hess_in)

% convert a tensor in matrix (second dim being constant or 1) 

for ii=1:size(hess_in,3)
    
   mat_out(:,ii) = hess_in(:,1,ii) ;
    
end