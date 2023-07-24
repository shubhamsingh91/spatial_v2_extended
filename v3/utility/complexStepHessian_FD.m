function [d2qdd_1, d2qdd_2,d2qdd_3] = complexStepHessian_FD(f, x, step)
    
% complex step for taking Hessian of Forward Dynamics

    if nargin == 2
        step = eps;
    end
    n = length(f(x));
    m = length(x);
    
    H = zeros(n,m,m);
    for ind = 1:m
       e = zeros(m,1);
       e(ind) = 1;

       [temp1, temp2 , temp3]= f(x+1i*e*step);
       
        d2qdd_1(:,:,ind) = imag( temp1)/step;  
        d2qdd_2(:,:,ind)= imag( temp2)/step;  
        d2qdd_3(:,:,ind)= imag( temp3)/step;  

    end
end