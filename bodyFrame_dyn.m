function  [bod] = bodyFrame_dyn( model, q, qd, qdd)
skew  = @(x) [ 0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0]; 
crm = @(x) [skew(x(1:3)), zeros(3,3) ; skew(x(4:6)) ,skew(x(1:3))]; 
a_grav = get_gravity(model);

if ~iscell(q) || ~iscell(qd) || ~iscell(qdd)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end

for i = 1:model.NB
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  vJ = S{i}*qd{i};
  Xup{i} = XJ * model.Xtree{i};
  if model.parent(i) == 0
    v{i} = vJ;
    a{i} = Xup{i}*(-a_grav) + S{i}*qdd{i};
  else
    v{i} = Xup{i}*v{model.parent(i)} + vJ;
    a{i} = Xup{i}*a{model.parent(i)} + S{i}*qdd{i} + crm(v{i})*vJ;
  end
  h{i} = model.I{i}*v{i};
  f{i} = model.I{i}*a{i} + crf(v{i})*h{i};
  
end

bod.f = f;

for i = model.NB:-1:1
  p = model.parent(i);
  if p ~= 0
        f{p} = f{p} + Xup{i}.'*f{i} ;
  end
end

bod.fC = f;
bod.v = f;
bod.a = a;

end