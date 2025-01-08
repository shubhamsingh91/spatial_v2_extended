function  [derivs] = spatial_force_SO_derivs(model, q, qd, qdd)
% SO Partial Derivatives for Cumulative spatial-force for multi-DoF joints
% Contributors - Shubham Singh, singh281@utexas.edu 
%
% This version caches:
%   - crf_bar(fCi) once per body i
%   - repeated crf(...) calls if used multiple times
%   - repeated IDOT(...) calls if used multiple times
% and preserves expr-# comments.

% -------------------------------------------------------------------------
% cross product matrix for force vectors with a bar over top
crf_bar = @(x)[ -skew(x(1:3)), -skew(x(4:6)); 
                -skew(x(4:6)),  zeros(3,3) ];

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd) || ~iscell(qdd)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end
if sum(model.has_rotor) > 1
    error('ID_derivatives does not support rotors');
end

% Gravity
a_grav = get_gravity(model);

% Copy inertia
IC = model.I;
I  = model.I;

% -------------------------------------------------------------------------
% Forward Pass
% -------------------------------------------------------------------------
for i = 1:model.NB
    [XJ, S{i}] = jcalc(model.jtype{i}, q{i});
    Xup{i}     = XJ * model.Xtree{i};
    
    if model.parent(i) == 0
        v{i}     = zeros(6,1);
        a{i}     = -a_grav;
        Xup0{i}  = Xup{i};
    else
        Xup0{i}  = Xup{i} * Xup0{model.parent(i)};
        v{i}     = v{model.parent(i)};
        a{i}     = a{model.parent(i)};
    end
    
    Xdown0{i} = inv(Xup0{i});
    S{i}      = Xdown0{i} * S{i};
    vJ{i}     = S{i} * qd{i};
    aJ{i}     = crm(v{i})*vJ{i} + S{i}*qdd{i};
    psid{i}   = crm(v{i}) * S{i};
    psidd{i}  = crm(a{i}) * S{i} + crm(v{i}) * psid{i};
    v{i}      = v{i} + vJ{i};
    a{i}      = a{i} + aJ{i};
    
    IC{i}     = Xup0{i}.' * I{i} * Xup0{i};
    Sd{i}     = crm(v{i}) * S{i};
    BC{i}     = 2 * factorFunctions(IC{i}, v{i});
    f{i}      = IC{i} * a{i} + crf(v{i}) * IC{i} * v{i};
    
    % Allocate partial derivatives
    d2fc_dq{i}  = zeros(6, model.NV, model.NV);
    d2fc_dv{i}  = zeros(6, model.NV, model.NV);
    d2fc_dav{i} = zeros(6, model.NV, model.NV);
    d2fc_dvq{i} = zeros(6, model.NV, model.NV);
end

% -------------------------------------------------------------------------
% Backward Pass
% -------------------------------------------------------------------------
for i = model.NB:-1:1

    BCi = BC{i};
    ICi = IC{i};
    fCi = f{i};

    % Cache crf_bar(fCi) once, since fCi is fixed for this body in the loop
    fCi_bar = crf_bar(fCi);

    for p = 1:model.nv(i)

        S_p     = S{i}(:,p);
        Sd_p    = Sd{i}(:,p);
        psid_p  = psid{i}(:,p);
        psidd_p = psidd{i}(:,p);

        % 2*I*C() patterns
        Bic_phii     = 2 * factorFunctions(ICi, S_p);
        Bic_psii_dot = 2 * factorFunctions(ICi, psid_p);
        
        % We'll see crf_bar(fCi)*S_p in multiple expressions, so cache it:
        fCi_Sp = fCi_bar * S_p;  % 6x1
        BCi_psid_p = BCi*psid_p;
        
        % Also cache IDOT(BCi,S_p), IDOT(ICi,S_p) if repeated:
        BCi_Sp = IDOT(BCi, S_p);  % 6x6
        ICi_Sp = IDOT(ICi, S_p);  % 6x6
        

        ii = model.vinds{i};
        j  = i;

        while j > 0
            jj = model.vinds{j};

            for t = 1:model.nv(j)

                S_t     = S{j}(:,t);
                Sd_t    = Sd{j}(:,t);
                psid_t  = psid{j}(:,t);
                psidd_t = psidd{j}(:,t);

                % 2*I*C() for S_t, psid_t
                Bic_phij      = 2 * factorFunctions(ICi, S_t);
                Bic_psijt_dot = 2 * factorFunctions(ICi, psid_t);
                
                % Also IDOT(BCi,S_t), IDOT(ICi,S_t)
                BCi_St = IDOT(BCi, S_t);
                ICi_St = IDOT(ICi, S_t);
                
                crfSt = crf(S_t);
                
                k = j;
                while k > 0
                    kk = model.vinds{k};

                    for r = 1:model.nv(k)

                        S_r     = S{k}(:,r);
                        Sd_r    = Sd{k}(:,r);
                        psid_r  = psid{k}(:,r);
                        psidd_r = psidd{k}(:,r);

                        % 2*I*C() for psid_r
                        Bic_psikr_dot = 2 * factorFunctions(ICi, psid_r);

                        % We do crf(S_r) multiple times:
                        crfSr = crf(S_r);
                        crmPsidr = crm(psid_r);

                        % ----------------------------------------------------------------
                        % expr-1 SO-q
                        % ----------------------------------------------------------------
                        % d2fc_dq{i}(:, jj(t), kk(r))
                        d2fc_dq{i}(:, jj(t), kk(r)) = ...
                             (Bic_psijt_dot + BCi_St)*psid_r ...
                             + ICi_St*psidd_r ...
                             + crfSr*(BCi*psid_t + ICi*psidd_t +  fCi_bar * S_t);

                        % ----------------------------------------------------------------
                        % expr-1 SO-av
                        % ----------------------------------------------------------------
                        d2fc_dav{i}(:, jj(t), kk(r)) = crfSr * (ICi*S_t);

                        % ----------------------------------------------------------------
                        % expr-1 SO-vq
                        % ----------------------------------------------------------------
                        d2fc_dvq{i}(:, jj(t), kk(r)) = ...
                            (Bic_psikr_dot + crfSr*BCi + 2*ICi*crmPsidr)*S_t ...
                            + crfSr*( ICi*(psid_t + Sd_t) );

                        % ================================================================
                        % if j ~= i => "off-diagonal" blocks
                        % ================================================================
                        if j ~= i
                            % expr-5 SO-q
                            % d2fc_dq{j}(:, kk(r), ii(p))
                            d2fc_dq{j}(:, kk(r), ii(p)) = ...
                                ICi_Sp*psidd_r ...
                                + crf_bar(ICi*psidd_p + fCi_Sp + BCi_psid_p)*S_r ...
                                + (Bic_psii_dot + BCi_Sp)*psid_r;

                            % expr-6 SO-q
                            d2fc_dq{j}(:, ii(p), kk(r)) = d2fc_dq{j}(:, kk(r), ii(p));

                            % expr-7 SO-v
                            d2fc_dv{j}(:, kk(r), ii(p)) = Bic_phii*S_r;
                            % expr-8 SO-v
                            d2fc_dv{j}(:, ii(p), kk(r)) = d2fc_dv{j}(:, kk(r), ii(p));

                            % expr-2 SO-av
                            d2fc_dav{j}(:, ii(p), kk(r)) = crfSr * (ICi*S_p);

                            % expr-5 SO-av
                            % d2fc_dav{j}(:, kk(r), ii(p))
                            d2fc_dav{j}(:, kk(r), ii(p)) = ICi_Sp * S_r;

                            % expr-2 SO-vq
                            d2fc_dvq{j}(:, ii(p), kk(r)) = ...
                                (Bic_psikr_dot + crfSr*BCi + 2*ICi*crmPsidr)*S_p ...
                                + crfSr*(ICi*(psid_p + Sd_p));

                            % expr-5 SO-vq
                            d2fc_dvq{j}(:, kk(r), ii(p)) = ...
                                (Bic_psii_dot + BCi_Sp)*S_r ...
                                + ICi_Sp*(psid_r + Sd_r);
                        end

                        % ================================================================
                        % if (k ~= j) => fill cross terms (kk < j)
                        % ================================================================
                        if (k ~= j)
                            % expr-2 SO-q
                            d2fc_dq{i}(:, kk(r), jj(t)) = d2fc_dq{i}(:, jj(t), kk(r));

                            % expr-3 SO-q
                            % d2fc_dq{k}(:, ii(p), jj(t))
                            d2fc_dq{k}(:, ii(p), jj(t)) = ...
                                  (Bic_psii_dot + BCi_Sp)*psid_t ...
                                  + ICi_Sp*psidd_t ...
                                  + crfSt*(BCi_psid_p + ICi*psidd_p + fCi_Sp);

                            % expr-1 SO-v
                            d2fc_dv{i}(:, jj(t), kk(r)) = Bic_phij*S_r;
                            % expr-2 SO-v
                            d2fc_dv{i}(:, kk(r), jj(t)) = d2fc_dv{i}(:, jj(t), kk(r));

                            % expr-3 SO-av
                            d2fc_dav{i}(:, kk(r), jj(t)) = ICi_St * S_r;

                            % expr-3 SO-vq
                            d2fc_dvq{i}(:, kk(r), jj(t)) = ...
                                (Bic_psijt_dot + BCi_St)*S_r ...
                                + ICi_St*(psid_r + Sd_r);

                            % ------------------------------------------------------------
                            % if (j ~= i) => (kk < j < i)
                            % ------------------------------------------------------------
                            if (j ~= i)
                                % expr-4 SO-q
                                d2fc_dq{k}(:, jj(t), ii(p)) = ...
                                    d2fc_dq{k}(:, ii(p), jj(t));

                                % expr-4 SO-v
                                d2fc_dv{k}(:, ii(p), jj(t)) = Bic_phii*S_t;
                                % expr-5 SO-v
                                d2fc_dv{k}(:, jj(t), ii(p)) = ...
                                    d2fc_dv{k}(:, ii(p), jj(t));

                                % expr-4 SO-av
                                d2fc_dav{k}(:, ii(p), jj(t)) = crfSt*ICi*S_p;

                                % expr-6 SO-av
                                d2fc_dav{k}(:, jj(t), ii(p)) = ICi_Sp * S_t;

                                % expr-4 SO-vq
                                tmp_c4 = BCi*S_p + ICi*(psid_p + Sd_p);
                                d2fc_dvq{k}(:, ii(p), jj(t)) = ...
                                    Bic_phii*psid_t + crf_bar(tmp_c4)*S_t;

                                % expr-6 SO-vq
                                d2fc_dvq{k}(:, jj(t), ii(p)) = ...
                                    (Bic_psii_dot + BCi_Sp)*S_t ...
                                    + ICi_Sp*(psid_t + Sd_t);

                            else
                                % (kk < j = i)
                                % expr-6 SO-v
                                d2fc_dv{k}(:, ii(p), jj(t)) = ...
                                    (crf(S_p)*ICi + crf_bar(ICi*S_p)) * S_t;
                            end

                        else
                            % (k == j => kk == jj(t))
                            % expr-3 SO-v
                            d2fc_dv{i}(:, jj(t), kk(r)) = ...
                                (crfSt*ICi + crf_bar(ICi*S_t)) * S_r;
                        end

                    end
                    k = model.parent(k);
                end
            end
            j = model.parent(j);
        end
    end

    % Accumulate to parent
    if model.parent(i) > 0
        p = model.parent(i);
        IC{p} = IC{p} + IC{i};
        BC{p} = BC{p} + BC{i};
        f{p}  = f{p}  + f{i};
    end
end

% Output partials
derivs.d2fc_dq  = d2fc_dq;
derivs.d2fc_dv  = d2fc_dv;
derivs.d2fc_dav = d2fc_dav;
derivs.d2fc_dvq = d2fc_dvq;

end % main function

%--------------------------------------------------------------------------
% IDOT(I, v) = crf(v)*I - I*crm(v)   (6x6)
%--------------------------------------------------------------------------
function IdotMat = IDOT(I, v)
    IdotMat = crf(v)*I - I*crm(v);
end
