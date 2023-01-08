function [derivs] = FD_SO_derivatives(model,q, qd,qdd)

[tau]      = ID(model, q ,qd ,qdd);                                        % Inverse dynamics
[dqdd_dq, dqdd_dv,~] = FD_derivatives( model, q, qd, tau );        % FD FO derivatives

[d2tau_dq, d2tau_dv, ddtau_dqv,M_FO] = ...
            ID_SO_derivatives_v7_pin_v5(model, q, qd,qdd);                 % ID SO derivatives
        
Hinv = Hinverse( model, q);                                                % Minv

derivs.d2FD_dq = mT(-Hinv,d2tau_dq+Tm(M_FO,dqdd_dq)+rotR(Tm(M_FO,dqdd_dq)));
derivs.d2FD_dv = mT(-Hinv,d2tau_dv);
derivs.d2FD_dqv = mT(-Hinv,ddtau_dqv+rotR(Tm(M_FO,dqdd_dv)));
derivs.d2FD_dtauq =  mT(-Hinv,Tm(M_FO,Hinv));

end 