function  derivs = modID_second_derivatives_soa( model, q, qd, qdd, lambda)
%%% Second derivatives of lambda^T tau, as needed in DDP

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if sum(model.has_rotor) > 1
    error('modID does not support rotors');
end

a_grav = get_gravity(model);
if ~iscell(q)
    [q, qd, qdd, lambda] = confVecToCell(model,q,qd,qdd, lambda);
end

Ic = model.I;

for i = 1:model.NB
  
  ii = model.vinds{i};
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
 
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
      X0{i} = Xup{i};
  else
      X0{i} = Xup{i}*X0{ model.parent(i) };
  end
  S{i} = X0{i}\S{i};
  vJ{i} = S{i}*qd{i};
  wJ{i} = S{i}*lambda{i};
 
  if model.parent(i) == 0
    a{i} = -a_grav;
    v{i} = zeros(6,1);
    w{i} = zeros(6,1);    
  else
    v{i} = v{model.parent(i)};
    w{i} = w{model.parent(i)};
    a{i} = a{model.parent(i)};
  end
  
  
  Yd{i}  = crm(v{i})*S{i};                         % = Sdot only when single-dof joints
  Ydd{i} = crm(a{i})*S{i}+crm(v{i})*Yd{i};         % Columns are of type spatial aDcelertion 
  Ud{i}  = 2*Yd{i}+crm(vJ{i})*S{i};                % Columns are of type spatial velocity
  Om{i}  = crm(w{i})*S{i};                         % = Yd when qdot = lambda
  
  v{i}  = v{i} + vJ{i};                            % Spatial velocithy
  a{i}  = a{i} + crm(v{i})*vJ{i} + S{i}*qdd{i};    % Spatial aDceleration
  w{i}  = w{i} + wJ{i};                            % Spatial velocity if qdot = lambda
   
  % Notes: Considering derviatives of vectors taken in local coordinates:
  % when j <= i (in a tree sense)
  % 1) dv_i / dq_j  = Yd{j}
  % 2) da_i / dq_j  = Ydd{j} - v{i} x Yd{j}
  % 3) dv_i / dqd_j = S{j} 
  % 4) da_i / dqd_j = Ud{j} - v{i} x S{j}
  
  Ic{i} = X0{i}'*Ic{i}*X0{i};                      % Initial composite inertia
  Bc{i} = 2*factorFunctions(Ic{i}, v{i});          % Initial composite coriolis term
  Dc{i} = 2*factorFunctions(Ic{i},-w{i});          % Unintuitive, just falls out of the math. Composite coriolis term when qdot = -lambda.
  
  h{i}  = Ic{i} *w{i};                             % Spatial momentum when qdot =lambda
  z{i}  = Bc{i}'*w{i};                             % Unintuitive, but has same "units" as spatial force
  f{i}  = Ic{i} *a{i} + crf(v{i})*Ic{i}*v{i};      % Spatial force
end

S    = myCell2Mat(model,S);
Yd   = myCell2Mat(model,Yd);
Ydd  = myCell2Mat(model,Ydd);
Ud   = myCell2Mat(model,Ud);
Om   = myCell2Mat(model,Om);

F1  = repmat(0*q{1}(1) , 6, model.NV);
F2  = repmat(0*q{1}(1) , 6, model.NV);
F3  = repmat(0*q{1}(1) , 6, model.NV);
F4  = repmat(0*q{1}(1) , 6, model.NV);
F5  = repmat(0*q{1}(1) , 6, model.NV);
F6  = repmat(0*q{1}(1) , 6, model.NV);
F7  = repmat(0*q{1}(1) , 6, model.NV);

dmod_dq  = repmat(0*q{1}(1), model.NV,1);
dmod_dv = repmat(0*q{1}(1), model.NV,1);
dtau_dq  = repmat(0*q{1}(1), model.NV,model.NV);
dtau_dv = repmat(0*q{1}(1), model.NV,model.NV);
dtau_da = repmat(0*q{1}(1), model.NV,model.NV);
dmod_dvv = repmat(0*q{1}(1), model.NV, model.NV);
dmod_dqq = repmat(0*q{1}(1), model.NV, model.NV);
dmod_dqv = repmat(0*q{1}(1), model.NV, model.NV);



for i = model.NB:-1:1
   ii = model.vinds{i};
   jj = model.subtree_vinds{i};
   kk = model.successor_vinds{i};
   
   % temp: only needed for this iteration of the loop
   hphi{i} = icrf(h{i})*S(:,ii);
   zphi{i} = icrf(z{i})*S(:,ii);
   fphi{i} = icrf(f{i})*S(:,ii);
   
   % Temp: used in subsequent iterations of the loop too
   F1(:,ii) = Ic{i} *S(:,ii);
   F2(:,ii) = Bc{i} *S(:,ii)  + Ic{i}*Ud(:,ii);
   F3(:,ii) = Bc{i} *Yd(:,ii) + Ic{i}*Ydd(:,ii) + fphi{i};
   F4(:,ii) = Bc{i}'*S(:,ii);
   F5(:,ii) = Dc{i}'*S(:,ii);
   F6(:,ii) = Dc{i} *Yd(:,ii) + Bc{i}'*Om(:,ii) + 2*icrf(h{i})*Yd(:,ii) + zphi{i};
   F7(:,ii) = Ic{i} *Om(:,ii) + hphi{i};
   
   % Notes: Considering derviatives of vectors taken in local coordinates:
   % when j > i (in a tree sense)
   % 1) df{i} / dq_j  = F3(:,j);   
   % 2) df{i} / dqd_j = F2(:,j);
   %
   % when j <= i (in a tree sense)
   % 3) df{i} / dq_j  = Bc{i}*Yd(:,j) + Ic{i}*Ydd(:,j); 
   % 4) df{i} / dqd_j = Ic{i}*Ud(:,j) + Bc{i}*S(:,j);
   %
   
   % :-)
   dmod_dq  (ii,1)   = Om(:,ii).'*f{i}  + Yd(:,ii).'  *z{i} + Ydd(:,ii).'*h{i};
   dmod_dv  (ii,1)   = S(:,ii).' *z{i}  + Ud(:,ii).'*h{i};   
   if length(kk) > 0
    dtau_dq  (ii,kk)  = S(:,ii).' *F3(:,kk);
    dtau_dv  (ii,kk)  = S(:,ii).' *F2(:,kk);
   end
   
   dtau_da  (ii,jj)  = S(:,ii).' *F1(:,jj);
   dmod_dvv (ii,jj)  = S(:,ii).' *F5(:,jj);
   dmod_dqv (ii,jj)  = Om(:,ii).'*F2(:,jj)  + Yd(:,ii).'*F5(:,jj);
   
   dtau_dq  (jj,ii)  = F1(:,jj).'*Ydd(:,ii) + F4(:,jj).'*Yd(:,ii);
   dtau_dv  (jj,ii)  = F1(:,jj).'*Ud(:,ii)  + F4(:,jj).'*S(:,ii);
   dmod_dqq (jj,ii)  = F3(:,jj).'*Om(:,ii)  + F6(:,jj).'*Yd(:,ii) + F7(:,jj).'*Ydd(:,ii);
   if length(kk) > 0
    dmod_dqv (kk,ii)  = F6(:,kk).'*S(:,ii)   + F7(:,kk).'*Ud(:,ii);
   end
   
   % Symmetry
   dmod_dvv (ii,ii)  = dmod_dvv(ii,ii)   + hphi{i}.' *S(:,ii); 
   if length(kk) > 0
    dmod_dvv (kk,ii)  = dmod_dvv(ii,kk).';
    dmod_dqq (ii,kk)  = dmod_dqq(kk,ii).';
    dtau_da  (kk,ii)  = dtau_da(ii,kk).';
   end
   
   p = model.parent(i);
   if p > 0
      z{p}  = z{p}  + z{i} ;
      h{p}  = h{p}  + h{i} ;
      f{p}  = f{p}  + f{i} ;
      Ic{p} = Ic{p} + Ic{i} ;
      Bc{p} = Bc{p} + Bc{i} ;
      Dc{p} = Dc{p} + Dc{i} ;
   end 
end

% return in struct
derivs.dtau_dq  = dtau_dq ;
derivs.dtau_dv  = dtau_dv;
derivs.dmod_dq  = dmod_dq ;
derivs.dmod_dv  = dmod_dv;
derivs.dmod_dqq = dmod_dqq;
derivs.dmod_dvv = dmod_dvv;
derivs.dmod_dqv = dmod_dqv ;

end

function M = myCell2Mat(model,S)
%     import casadi.*
%     M = SX.zeros(6,model.NV);
%     M = zeros(6,model.NV);
    for i = 1:length(S)
        M(:,model.vinds{i}) = S{i};
    end
end