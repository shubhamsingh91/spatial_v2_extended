% SecondOrderDerivation.m
% -----------------------------------------------------------------------------
% Complex-step verification of the global-notation second-order BUILDING-BLOCK
% identities from the appendix of the SO rigid-body-dynamics derivatives paper.
%
% For a fixed (stacked, per-body) force vector Y, with  S_Y = Y (x)* S  and
% Psid_Y = Y (x)* Psid  (where (x)* is the icrf operator applied block-wise):
%
%   grad_q   Psid'  Y = S_Y' * Sig * Psid
%   grad_q   Psidd' Y = S_Y' * Sig * Psidd + 2 * Psid_Y' * Sig * Psid
%   grad_q   Om'    Y = S_Y' * Sig * Om
%   grad_qd  Psid'  Y = S_Y' * Sigbar * S
%   grad_qd  Upd'   Y = S_Y' * (Sig + Sigbar) * S
%   grad_qd  Psidd' Y = S_Y' * Sigbar * Upd + 2 * Psid_Y' * Sigbar * S
%
% Global notation (body coordinates):
%   S      = BlkDiag(S_1,...,S_N)                 (6N x NV, body i in its block-row)
%   Sig    = predecessor-sum operator, Sig(i,j) = ^iX_j for j <= i, I on diagonal
%   Sigbar = Sig - I
%   Psid_i = Vp_i x S_i,  Psidd_i = Ap_i x S_i + Vp_i x (Vp_i x S_i)
%   Upd_i  = v_i x S_i + Vp_i x S_i,   Om_i = Wp_i x S_i   (w-recursion uses lambda)
%
% Tested over many random trees / sizes for robustness.
% -----------------------------------------------------------------------------

clear; clc;
NTRIALS = 25;
TOL     = 1e-7;
names   = {'grad_q  Psid','grad_q  Psidd','grad_q  Om', ...
           'grad_qd Psid','grad_qd Upd','grad_qd Psidd'};
maxerr  = zeros(1,6);

fprintf('=====================================================\n');
fprintf(' Global-notation SO building-block identities (cstep)\n');
fprintf(' %d random trials, sizes N = 1..8\n', NTRIALS);
fprintf('=====================================================\n');

for t = 1:NTRIALS
    rng(t);                              % reproducible per trial
    N = randi([1 8]);
    model = autoTree(N, 1.5, pi/3);
    model = postProcessModel(model);
    for i = 1:model.NB
        model.I{i} = inertiaVecToMat( rand(10,1) );
    end
    e = checkIdentities(model);
    maxerr = max(maxerr, e(:).');
end

for k = 1:6
    fprintf('%-14s  max residual = %.3e\n', names{k}, maxerr(k));
end
fprintf('-----------------------------------------------------\n');
if max(maxerr) > TOL
    error('A global-notation identity is out of tolerance (max %.3e)', max(maxerr));
end
fprintf(' ALL %d TRIALS PASSED  (overall max residual = %.2e)\n', NTRIALS, max(maxerr));
fprintf('=====================================================\n');


% ============================ verification ==========================

function e = checkIdentities(model)
    q      = normalizeConfVec(model, rand(model.NQ,1));
    qd     = rand(model.NV,1);
    qdd    = rand(model.NV,1);
    lambda = rand(model.NV,1);
    Y      = rand(6*model.NB,1);     % fixed stacked force vector

    G            = buildGlobals(model, q, qd, qdd, lambda);
    [Sig,Sigbar] = buildSigma(model, q);
    SY           = blockIcrf(model, Y, G.S);       % S_Y    = Y (x)* S
    PsidY        = blockIcrf(model, Y, G.Psid);    % Psid_Y = Y (x)* Psid

    % grad_q identities (perturb q)
    gq_Psid  = complexStepJacobian(@(x) tY(gfield(model,x,qd,qdd,lambda,'Psid' ), Y), q);
    gq_Psidd = complexStepJacobian(@(x) tY(gfield(model,x,qd,qdd,lambda,'Psidd'), Y), q);
    gq_Om    = complexStepJacobian(@(x) tY(gfield(model,x,qd,qdd,lambda,'Om'   ), Y), q);

    % grad_qd identities (perturb qd)
    gqd_Psid  = complexStepJacobian(@(x) tY(gfield(model,q,x,qdd,lambda,'Psid' ), Y), qd);
    gqd_Upd   = complexStepJacobian(@(x) tY(gfield(model,q,x,qdd,lambda,'Upd'  ), Y), qd);
    gqd_Psidd = complexStepJacobian(@(x) tY(gfield(model,q,x,qdd,lambda,'Psidd'), Y), qd);

    e = [ rerr(gq_Psid , SY.'*Sig*G.Psid)
          rerr(gq_Psidd, SY.'*Sig*G.Psidd + 2*PsidY.'*Sig*G.Psid)
          rerr(gq_Om   , SY.'*Sig*G.Om)
          rerr(gqd_Psid , SY.'*Sigbar*G.S)
          rerr(gqd_Upd  , SY.'*(Sig+Sigbar)*G.S)
          rerr(gqd_Psidd, SY.'*Sigbar*G.Upd + 2*PsidY.'*Sigbar*G.S) ];
end


% ============================ helpers ===============================

function out = tY(M, Y)
    out = M.' * Y;          % out(k) = M(:,k).' * Y   (non-conjugate transpose)
end

function v = rerr(A, B)
    v = norm(A(:) - B(:));
end

function M = gfield(model, q, qd, qdd, lambda, fname)
    G = buildGlobals(model, q, qd, qdd, lambda);
    M = G.(fname);
end

function G = buildGlobals(model, q, qd, qdd, lambda)
    if ~iscell(q)
        [q, qd, qdd, lambda] = confVecToCell(model, q, qd, qdd, lambda);
    end
    a_grav = get_gravity(model);
    NB = model.NB; NV = model.NV;
    S = zeros(6*NB, NV); Psid = S; Psidd = S; Upd = S; Om = S;
    v = cell(NB,1); a = cell(NB,1); w = cell(NB,1);
    for i = 1:NB
        ii = model.vinds{i};
        ri = 6*(i-1) + (1:6);
        [XJ, Si] = jcalc(model.jtype{i}, q{i});
        Xup = XJ * model.Xtree{i};
        vJ = Si*qd{i};
        wJ = Si*lambda{i};
        if model.parent(i) == 0
            Vp = zeros(6,1); Ap = Xup*(-a_grav); Wp = zeros(6,1);
        else
            Vp = Xup*v{model.parent(i)};
            Ap = Xup*a{model.parent(i)};
            Wp = Xup*w{model.parent(i)};
        end
        v{i} = Vp + vJ;
        a{i} = Ap + crm(v{i})*vJ + Si*qdd{i};
        w{i} = Wp + wJ;
        S(ri,ii)     = Si;
        Psid(ri,ii)  = crm(Vp)*Si;
        Psidd(ri,ii) = crm(Ap)*Si + crm(Vp)*(crm(Vp)*Si);
        Upd(ri,ii)   = crm(v{i})*Si + crm(Vp)*Si;
        Om(ri,ii)    = crm(Wp)*Si;
    end
    G.S = S; G.Psid = Psid; G.Psidd = Psidd; G.Upd = Upd; G.Om = Om;
end

function [Sig, Sigbar] = buildSigma(model, q)
    if ~iscell(q)
        z = zeros(model.NV,1);
        [q, ~, ~, ~] = confVecToCell(model, q, z, z, z);
    end
    NB = model.NB;
    X0 = cell(NB,1);
    for i = 1:NB
        [XJ,~] = jcalc(model.jtype{i}, q{i});
        Xup = XJ * model.Xtree{i};
        if model.parent(i) == 0, X0{i} = Xup; else, X0{i} = Xup*X0{model.parent(i)}; end
    end
    Sig = zeros(6*NB);
    for i = 1:NB
        j = i;
        while j > 0
            Sig(6*(i-1)+(1:6), 6*(j-1)+(1:6)) = X0{i}/X0{j};   % ^iX_j
            j = model.parent(j);
        end
    end
    Sigbar = Sig - eye(6*NB);
end

function MY = blockIcrf(model, Y, M)
    % MY = Y (x)* M : block-diagonal application of icrf(Y_i) to body i's block
    NB = model.NB;
    MY = zeros(size(M));
    for i = 1:NB
        ii = model.vinds{i};
        ri = 6*(i-1) + (1:6);
        MY(ri,ii) = icrf(Y(ri)) * M(ri,ii);
    end
end
