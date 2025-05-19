function derivs = modID_second_derivatives_ground(model,q,qd,qdd,lambda)
% Forward-mode accumulation of second derivatives for the spatial
% articulated-body inverse-dynamics model (ground contact version).
% The per-link 6×NV derivative blocks are flattened into one matrix
% of size 6×(NV*NB).  Otherwise identical to the original code.

% -------------------------------------------------------------------------
% bookkeeping
% -------------------------------------------------------------------------
if ~isfield(model,'nq'),   model = postProcessModel(model);  end
if sum(model.has_rotor) > 1
    error('modID does not support rotors');
end

NV  = model.NV;            % total DoF
NB  = model.NB;            % number of bodies
col = @(i) ((i-1)*NV + 1):(i*NV);   % columns belonging to body i

% -------------------------------------------------------------------------
% gravity & configuration reshaping
% -------------------------------------------------------------------------
a_grav = get_gravity(model);
if ~iscell(q)
    [q,qd,qdd,lambda] = confVecToCell(model,q,qd,qdd,lambda);
end

% -------------------------------------------------------------------------
% flat derivative blocks (6 × NV·NB)
% -------------------------------------------------------------------------
dv_dq_p  = zeros(6,NV*NB);   dv_dqd_p = zeros(6,NV*NB);
da_dq_p  = zeros(6,NV*NB);   dw_dq_p  = zeros(6,NV*NB);

dv_dq    = zeros(6,NV*NB);   dv_dqd   = zeros(6,NV*NB);
da_dq    = zeros(6,NV*NB);   dw_dq    = zeros(6,NV*NB);

dh_dq    = zeros(6,NV*NB);   dz_dq    = zeros(6,NV*NB);
dz_dqd   = zeros(6,NV*NB);   df_dq    = zeros(6,NV*NB);

% Hessians
H_qdq  = zeros(NV,NV);
H_qdqd = zeros(NV,NV);
H_qq   = zeros(NV,NV);

% -------------------------------------------------------------------------
% per-link spatial variables (unchanged cell layout)
% -------------------------------------------------------------------------
S  = cell(NB,1);   Xup = cell(NB,1);
v  = cell(NB,1);   a   = cell(NB,1);   w  = cell(NB,1);
vp = cell(NB,1);   ap  = cell(NB,1);   wp = cell(NB,1);

h  = cell(NB,1);   z   = cell(NB,1);   f  = cell(NB,1);

I = model.I;

% =========================================================================
% 1) forward (outward) pass
% =========================================================================
for i = 1:NB
    ii    = model.vinds{i};        % global indices of this joint
    idx_i = col(i);                % columns for body i
    c_ii  = (i-1)*NV + ii;         % joint-specific slice

    % joint kinematics
    [XJ,S{i}] = jcalc(model.jtype{i},q{i});
    vJ  = S{i}*qd{i};
    wJ  = S{i}*lambda{i};
    Xup{i} = XJ * model.Xtree{i};

    if model.parent(i)==0                                     % ROOT
        vp{i} = zeros(6,1);           wp{i} = zeros(6,1);
        ap{i} = Xup{i} * (-a_grav);

        da_dq_p(:,c_ii) = crm(ap{i}) * S{i};
    else                                                       % CHILD
        p      = model.parent(i);
        idx_p  = col(p);

        vp{i}  = Xup{i} * v{p};
        wp{i}  = Xup{i} * w{p};
        ap{i}  = Xup{i} * a{p};

        % propagate derivatives from parent
        dv_dq_p(:,idx_i)   = Xup{i} * dv_dq(:,idx_p);
        da_dq_p(:,idx_i)   = Xup{i} * da_dq(:,idx_p);
        dw_dq_p(:,idx_i)   = Xup{i} * dw_dq(:,idx_p);
        dv_dqd_p(:,idx_i)  = Xup{i} * dv_dqd(:,idx_p);

        % plus local contributions
        dv_dq_p(:,c_ii) = crm(vp{i}) * S{i};
        da_dq_p(:,c_ii) = crm(ap{i}) * S{i};
        dw_dq_p(:,c_ii) = crm(wp{i}) * S{i};
    end

    % spatial velocities / accelerations
    v{i} = vp{i} + vJ;
    a{i} = ap{i} + crm(v{i})*vJ + S{i}*qdd{i};
    w{i} = wp{i} + wJ;

    dv_dq(:,idx_i)  = dv_dq_p(:,idx_i);

    da_dq(:,idx_i)  = da_dq_p(:,idx_i) - crm(vJ) * dv_dq(:,idx_i);

    dw_dq(:,idx_i)  = dw_dq_p(:,idx_i);

    dv_dqd(:,idx_i) = dv_dqd_p(:,idx_i);
    dv_dqd(:,c_ii)  = S{i};

    % momentum-like terms
    h{i}            = I{i} * w{i};
    dh_dq(:,idx_i)  = I{i} * dw_dq(:,idx_i);

    z{i}            = I{i}*crm(w{i})*v{i} - crf(w{i})*I{i}*v{i};

    dz_dq(:,idx_i)  = I{i}*crm(w{i})*dv_dq(:,idx_i) ...
                    - I{i}*crm(v{i})*dw_dq(:,idx_i) ...
                    - crf(w{i})*I{i}*dv_dq(:,idx_i) ...
                    - icrf(I{i}*v{i})*dw_dq(:,idx_i);

    dz_dqd(:,idx_i) = I{i}*crm(w{i})*dv_dqd(:,idx_i) ...
                    - crf(w{i})*I{i}*dv_dqd(:,idx_i);

    f{i}            = I{i}*a{i} + crf(v{i})*I{i}*v{i};
    df_dq(:,idx_i)  = I{i}*da_dq(:,idx_i) ...
                    + crf(v{i})*I{i}*dv_dq(:,idx_i) ...
                    + icrf(I{i}*v{i})*dv_dq(:,idx_i);
end

% =========================================================================
% 2) backward (inward) pass
% =========================================================================
for i = NB:-1:1
    ii    = model.vinds{i};
    idx_i = col(i);      c_ii = (i-1)*NV + ii;

    z{i}           = z{i} + crf(S{i}*qd{i}) * h{i};

    dz_dq(:,idx_i) = dz_dq(:,idx_i) + crf(S{i}*qd{i}) * dh_dq(:,idx_i);
    dz_dqd(:,c_ii) = dz_dqd(:,c_ii) + icrf(h{i}) * S{i};

    % Hessian rows
    H_qdq(ii,:)  = S{i}' * ( dz_dq(:,idx_i) ...
                  - crf(v{i}) * dh_dq(:,idx_i) ...
                  - icrf(h{i}) * dv_dq(:,idx_i) );

    H_qdqd(ii,:) = S{i}' * ( dz_dqd(:,idx_i) ...
                  - icrf(h{i}) * dv_dqd(:,idx_i) );

    H_qq(ii,:)   = -S{i}' * ( crf(vp{i}) * dz_dq(:,idx_i) ...
                  + crf(ap{i}) * dh_dq(:,idx_i) ...
                  + crf(wp{i}) * df_dq(:,idx_i) ...
                  + icrf(z{i}) * dv_dq_p(:,idx_i) ...
                  + icrf(h{i}) * da_dq_p(:,idx_i) ...
                  + icrf(f{i}) * dw_dq_p(:,idx_i) );

    % propagate to parent
    p = model.parent(i);
    if p > 0
        idx_p = col(p);

        z{p}           = z{p} + Xup{i}' * z{i};

        dz_dq(:,idx_p) = dz_dq(:,idx_p) + Xup{i}' * dz_dq(:,idx_i);
        dz_dq(:, (p-1)*NV + ii) = dz_dq(:, (p-1)*NV + ii) ...
                                + Xup{i}' * icrf(z{i}) * S{i};

        dz_dqd(:,idx_p) = dz_dqd(:,idx_p) + Xup{i}' * dz_dqd(:,idx_i);

        h{p}           = h{p} + Xup{i}' * h{i};

        dh_dq(:,idx_p) = dh_dq(:,idx_p) + Xup{i}' * dh_dq(:,idx_i);
        dh_dq(:, (p-1)*NV + ii) = dh_dq(:, (p-1)*NV + ii) ...
                                + Xup{i}' * icrf(h{i}) * S{i};

        f{p}           = f{p} + Xup{i}' * f{i};

        df_dq(:,idx_p) = df_dq(:,idx_p) + Xup{i}' * df_dq(:,idx_i);
        df_dq(:, (p-1)*NV + ii) = df_dq(:, (p-1)*NV + ii) ...
                                + Xup{i}' * icrf(f{i}) * S{i};
    end
end

% -------------------------------------------------------------------------
% results
% -------------------------------------------------------------------------
derivs.dmod_dqq = H_qq;
derivs.dmod_dvv = H_qdqd;
derivs.dmod_dqv = H_qdq.';
end
