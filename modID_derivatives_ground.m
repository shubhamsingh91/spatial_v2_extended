function  [grad_q, grad_qd, grad_qdd, total_grad] = modID_derivatives_ground( model, q, qd, qdd, lambda )

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
a_grav = get_gravity(model);
IC = model.I;
I = model.I;

% Calculate the derivatives using reverse mode chain rule
for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  Xup{i} = XJ * model.Xtree{i}; % i_X_p(i)
  
  if model.parent(i) == 0
    vp{i} = zeros(6,1);
    ap{i} = -a_grav;
    wp{i} = zeros(6,1);
    Xup0{i} = Xup{i}; % i_X_0
  else
    Xup0{i} = Xup{i}*Xup0{model.parent(i)}; % i_X_0
    vp{i} = v{model.parent(i)};
    wp{i} = w{model.parent(i)};
    ap{i} = a{model.parent(i)};
  end
  Xdown0{i} = inv(Xup0{i}); %0_X_i
  
   S{i} = Xdown0{i}*S{i}; %0_S_i
   vJ{i} = S{i}*qd{i};
   wJ{i} = S{i}*lambda{i};
  
  Yd{i}   = crm(vp{i})*S{i};
  Ydd{i}  = crm(vp{i})*Yd{i} + crm(ap{i})*S{i};
  Psid{i} = 2*Yd{i}+crm(vJ{i})*S{i};

  
  v{i} = vp{i} + vJ{i};
  a{i} = ap{i} + crm(v{i})*vJ{i} + S{i}*qdd{i};
  IC{i} = Xup0{i}.'*I{i}*Xup0{i};

  f{i} = IC{i}*a{i} + crf(v{i})*IC{i}*v{i};
  
  w{i} = wp{i} + wJ{i};
  h{i} = IC{i}*w{i}; 
  B{i} = factorFunctions(IC{i},v{i});
  z{i} = B{i}.'*w{i};
end

grad_q = repmat( 0*q{1}(1), [size(qd,1),size(lambda,2)] );
grad_qd = repmat(0*0*q{1}(1),[size(qd,1),size(lambda,2)] );
grad_qdd = repmat(0*0*q{1}(1),[size(qdd,1),size(lambda,2)] );

for i = model.NB:-1:1
   ii = model.vinds{i};    
   grad_qd(ii,:) = Psid{i}.'*h{i} + 2*S{i}.'*z{i};
   grad_q(ii,:)  = (icrf(f{i})*S{i}).'*wp{i} + 2*Yd{i}.'*z{i} + Ydd{i}.'*h{i};
   grad_qdd(ii,:) = S{i}.'*h{i};

   p = model.parent(i);
   if p > 0
      z{p} = z{p}+z{i};
      h{p} = h{p}+h{i};
      f{p} = f{p}+f{i};
   end 
end
total_grad = [ grad_q ; grad_qd; grad_qdd];
