function [derivs] = FD_SO_derivatives(model,q, qd,qdd)

[tau]      = ID(model, q ,qd ,qdd);                                        % Inverse dynamics
[dqdd_dq, dqdd_dv,~] = FD_derivatives( model, q, qd, tau );        % FD FO derivatives

derivs_IDSO =  ID_SO_derivatives(model, q, qd,qdd);                   % ID SO derivatives
        
d2tau_dq=derivs_IDSO.d2tau_dq;
d2tau_dv=derivs_IDSO.d2tau_dqd;
ddtau_dqv=derivs_IDSO.d2tau_cross;
M_FO=derivs_IDSO.dM_dq;


Hinv = Hinverse( model, q);                                                % Minv

derivs.d2FD_dq = mT(-Hinv,d2tau_dq+Tm(M_FO,dqdd_dq)+rotR(Tm(M_FO,dqdd_dq)));
derivs.d2FD_dv = mT(-Hinv,d2tau_dv);
derivs.d2FD_dqv = mT(-Hinv,ddtau_dqv+rotR(Tm(M_FO,dqdd_dv)));
derivs.d2FD_dtauq =  mT(-Hinv,Tm(M_FO,Hinv));

end 