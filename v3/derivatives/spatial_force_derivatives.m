function  [df_dq,dfc_dq,df_dv,dfc_dv,df_da,dfc_da] = spatial_force_derivatives( model, q, qd, qdd, i_idx,j_idx )

if ~isfield(model,'nq')
    model = postProcessModel(model);
end
if ~iscell(q) || ~iscell(qd) || ~iscell(qdd)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end
if sum(model.has_rotor) > 1
    error('ID_derivatives does not support rotors');
end

a_grav = get_gravity(model);
IC = model.I;
I = model.I;
I_0 = model.I;

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  Xup{i} = XJ * model.Xtree{i};
  
  if model.parent(i) == 0
    v{i}  = zeros(6,1);
    a{i}  = -a_grav;
    Xup0{i} = Xup{i};
  else
    Xup0{i} = Xup{i}*Xup0{model.parent(i)};
    v{i}  = v{model.parent(i)};
    a{i}  = a{model.parent(i)};
  end
  
  Xdown0{i} = inv(Xup0{i});
  
  S{i} = Xdown0{i}*S{i};
  vJ{i}= S{i}*qd{i};
  aJ{i} = crm(v{i})*vJ{i} + S{i}*qdd{i};
  
  psid{i} = crm(v{i})*S{i};
  psidd{i}= crm(a{i})*S{i} + crm(v{i})*psid{i};

  v{i} = v{i} + vJ{i};
  a{i} = a{i} + aJ{i};
  IC{i} = Xup0{i}.'*I{i}*Xup0{i};
  I_0{i} = Xup0{i}.'*I_0{i}*Xup0{i};
  Sd{i} = crm(v{i})*S{i};

  BC{i} = factorFunctions(IC{i},v{i});
  B{i} = factorFunctions(I_0{i},v{i});

  f{i}  =  IC{i}*a{i} + crf(v{i})*IC{i}*v{i};
  fc{i} = f{i};
  
end


for i = model.NB:-1:1

  if model.parent(i) > 0
     p = model.parent(i);
     IC{p} = IC{p} + IC{i};
     BC{p} = BC{p} + BC{i};
     fc{p}  = fc{p}  + fc{i};
  end
  
end

jj = model.vinds{j_idx}; ii = model.vinds{i_idx};

df_dq = zeros(6,model.nv(j_idx));
dfc_dq = zeros(6,model.nv(j_idx));

df_dv = zeros(6,model.nv(j_idx));
dfc_dv = zeros(6,model.nv(j_idx));

df_da = zeros(6,model.nv(j_idx));
dfc_da = zeros(6,model.nv(j_idx));

if ismember(j_idx,model.ancestors{i_idx}) 
    
    df_dq = I_0{i_idx}*psidd{j_idx} + cmf_bar(f{i_idx})*S{j_idx} + 2*B{i_idx}*psid{j_idx};
    dfc_dq = IC{i_idx}*psidd{j_idx} + cmf_bar(fc{i_idx})*S{j_idx} + 2*BC{i_idx}*psid{j_idx};
    
    df_dv = I_0{i_idx}*(Sd{j_idx} + psid{j_idx}) + 2*B{i_idx}*S{j_idx};
    dfc_dv = IC{i_idx}*(Sd{j_idx} + psid{j_idx}) + 2*BC{i_idx}*S{j_idx};
    
    df_da = I_0{i_idx}*S{j_idx};
    dfc_da = IC{i_idx}*S{j_idx};

elseif (ismember(i_idx,model.ancestors{j_idx}))
        
    dfc_dq =  IC{j_idx}*psidd{j_idx} + cmf_bar(fc{j_idx})*S{j_idx} + 2*BC{j_idx}*psid{j_idx};
    dfc_dv = IC{j_idx}*(Sd{j_idx} + psid{j_idx}) + 2*BC{j_idx}*S{j_idx};
    
    dfc_da = IC{j_idx}*S{j_idx};

end



    
end
