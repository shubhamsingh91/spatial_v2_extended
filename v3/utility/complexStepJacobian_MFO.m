function [H, evals, steps] = complexStepJacobian_MFO(f, x, step)

   
% This calculates the jacobian of M matrix (or M FO)

    if nargin == 2
        step = eps;
    end
    
    n = size(f(x),1);
    m = length(x);
    
    H = zeros(n,m,m);
    
    for ind = 1:m
       e = zeros(m,1);
       e(ind) = 1;
       steps(:,ind) = x+1i*e*step;
       evals(:,:,ind) = f(x+1i*e*step);
       
       H(:,:,ind) = imag( f(x+1i*e*step))/step;  
    end
end