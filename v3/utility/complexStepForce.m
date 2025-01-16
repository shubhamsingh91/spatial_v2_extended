function [d2fc_dq_g,d2fc_dq_b] = complexStepForce(model, f, x, jj,kk, step)
    if nargin == 5
        step = eps;
    end
    n = 6;    
    d2fc_dq_g = zeros(6,model.nv(jj),model.nv(kk));
    d2fc_dq_b = zeros(6,model.nv(jj),model.nv(kk));

    m = length(x);
    
    k_idx_begin = model.vinds{kk}(1);
    k_idx_end = model.vinds{kk}(end);

    for i = 1:model.nv(kk)
        
       e = zeros(m,1);
       e(i+k_idx_begin-1) = 1;
       [df_dq,dfci_dqj,df_dv,dfc_dv,df_da,dfc_da,dfci_dqj_body,dfc_dv_body,dfc_da_body]= f(x+1i*e*step);
       
        d2fc_dq_g(:,:,i) = imag( dfci_dqj)/step;  % d2fci_dqj_dqk
        d2fc_dq_b(:,:,i) = imag( dfci_dqj_body)/step;  % d2fci_dqj_dqk

    end
end