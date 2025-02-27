function [d2fc_dadq,d2fc_dadq_b] = complexStepForce_SOaq(model, f, x, jj,kk, step)
    if nargin == 5
        step = eps;
    end
    n = 6;    
    d2fc_dadq = zeros(6,model.nv(jj),model.nv(kk));
    d2fc_dadq_b = zeros(6,model.nv(jj),model.nv(kk));

    m = length(x);
    
    k_idx_begin = model.vinds{kk}(1);
    k_idx_end = model.vinds{kk}(end);

    for i = 1:model.nv(kk)
        
       e = zeros(m,1);
       e(i+k_idx_begin-1) = 1;
       [df_dq,dfci_dqj,df_dv,dfci_dvj,df_da,dfci_daj,dfci_dqj_b,dfci_dvj_b,dfci_daj_b]= f(x+1i*e*step); % returns df/da
       
        d2fc_dadq(:,:,i) = imag( dfci_daj)/step;  % d/dq(df/da) = d2fi_daj_dqk
        d2fc_dadq_b(:,:,i) = imag( dfci_daj_b)/step;  % d/dq(df/da) = d2fi_daj_dqk

    end
end