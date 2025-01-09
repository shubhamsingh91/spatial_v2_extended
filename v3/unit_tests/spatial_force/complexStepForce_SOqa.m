function [d2fc_dqda] = complexStepForce_SOqa(model, f, x, jj,kk, step)
    if nargin == 5
        step = eps;
    end
    n = 6;    
    d2fc_dqda = zeros(6,model.nv(jj),model.nv(kk));
    m = length(x);
    
    k_idx_begin = model.vinds{kk}(1);
    k_idx_end = model.vinds{kk}(end);

    for i = 1:model.nv(kk)
        
       e = zeros(m,1);
       e(i+k_idx_begin-1) = 1;
       [df_dq,dfci_dqj,df_dv,dfci_dvj,df_da,dfc_da]= f(x+1i*e*step); % returns df/da
       
        d2fc_dqda(:,:,i) = imag( dfci_dqj)/step;  % d/da(df/dq) = d2fi_dqj_dak

    end
end