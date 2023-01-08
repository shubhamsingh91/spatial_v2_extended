function out = tens_norm(hess_in)

for ii=1:size(hess_in,3)
    test1(ii) = norm(hess_in(:,:,ii));
end

out = norm(test1);