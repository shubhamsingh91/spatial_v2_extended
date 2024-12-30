function [d2fc_dadq] = complexStepForce_SOaq(model, f, x, jj,kk, step)
    if nargin == 5
        step = eps;
    end
    n = 6;    
    d2fc_dadq = zeros(6,model.nv(jj),model.nv(kk));
    m = length(x);
    
    k_idx_begin = model.vinds{kk}(1);
    k_idx_end = model.vinds{kk}(end);

    for i = 1:model.nv(kk)
        
       e = zeros(m,1);
       e(i+k_idx_begin-1) = 1;
       [df_dq,dfci_dqj,df_dv,dfci_dvj,df_da,dfci_daj]= f(x+1i*e*step); % returns df/da
       
        d2fc_dadq(:,:,i) = imag( dfci_daj)/step;  % d/dq(df/da) = d2fi_daj_dqk

    end
end