function  [glob, bod] = GlobalDynamics_red( model, q, qd, qdd)
skew  = @(x) [ 0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0]; 
crm = @(x) [skew(x(1:3)), zeros(3,3) ; skew(x(4:6)) ,skew(x(1:3))]; 

if ~iscell(q) || ~iscell(qd) || ~iscell(qdd)
    [q, qd, qdd] = confVecToCell(model,q,qd,qdd);
end

a_grav = get_gravity(model);
IC = model.I;
I = model.I;

for i = 1:model.NB
    
  [ XJ, S{i} ] = jcalc( model.jtype{i}, q{i} );
  Xup{i} = XJ * model.Xtree{i};
  
     if model.parent(i)==0
            X_m_0{i} = Xup{i}; 
            X_0_m{i} = inv(Xup{i});
      else
            X_m_0{i} = Xup{i}*X_m_0{model.parent(i)}; 
            X_0_m{i} = X_0_m{model.parent(i)}*inv(Xup{i});
     end 
     
     IC{i} = X_m_0{i}.'*model.I{i}*X_m_0{i}; % IC in ground frame
     I{i} =  X_m_0{i}.'*model.I{i}*X_m_0{i};  % I in ground frame
     
     S{i} = inv(X_m_0{i})*S{i};                 % phi in ground frame
     vJ{i} = S{i}*qd{i};
      
  if model.parent(i) == 0
    v{i}  = vJ{i};
    a{i}  = (-a_grav) + crm(v{i})*vJ{i} + S{i}*qdd{i};
  else
    v{i}  = v{model.parent(i)} + vJ{i};
    a{i}  = a{model.parent(i)} + crm(v{i})*vJ{i} + S{i}*qdd{i};
  end
  aJ = crm(v{i})*vJ{i} + S{i}*qdd{i}; 
  % spatial kinematic quantities in ground frame
  
  Sd{i} = crm(v{i})*S{i};                                   % phi_dot
  psid{i} = crm(v{i}-vJ{i})*S{i};                           % psi_dot
  psidd{i}= crm(a{i}-aJ)*S{i} + crm(v{i}-vJ{i})* psid{i};   % psi_ddot
   
  BC{i} = factorFunctions(I{i},v{i});
  f{i}  =  I{i}*a{i} + crf(v{i})*I{i}*v{i};
  vJxS{i} = crm(vJ{i})*S{i}; % vJi x S{i}
  Ia{i} = I{i}*a{i};
  
  % Body coordinates
  Sb{i} =   X_m_0{i}*S{i};
  psidb{i} =X_m_0{i}*psid{i};
  Sdb{i} =  X_m_0{i}*Sd{i};
  vb{i} =   X_m_0{i}*v{i};
  vJxSb{i}=X_m_0{i}*vJxS{i};
  psiddb{i} =X_m_0{i}*psidd{i};
  Ib{i} =   model.I{i};
  ab{i} =   X_m_0{i}*a{i};
  Iab{i} =   X_m_0{i}*Ia{i};
  fb{i} =   X_m_0{i}*f{i};  % is this the correct transform?

end

glob.B = BC;  % B(I,v) 
glob.f = f;  % fi

for i = model.NB:-1:1
    
  if model.parent(i) > 0
     p = model.parent(i);
     IC{p} = IC{p} + IC{i};
     BC{p} = BC{p} + BC{i};
     f{p}  = f{p}  + f{i};

  end 
    
end

for i=1:model.NB
    
    ICb{i} =  X_m_0{i}*IC{i}* X_0_m{i};
       fCb{i} =  X_0_m{i}*f{i};
 
end


  % glob coordinates
  glob.S = S;           glob.Sd = Sd;      glob.psid = psid;
  glob.psidd = psidd;   glob.BC = BC;      glob.fC = f;
  glob.v = v;           glob.a = a;        glob.vJ = vJ;  
  glob.vJxS = vJxS;     glob.I = I;        glob.IC = IC;
  glob.Ia = Ia;
  glob.iXO = X_m_0; 
  glob.OXi = X_0_m;
  
  % body coordinates
  bod.S = Sb;  bod.Sd = Sdb;    bod.psid = psidb;
  bod.v = vb;  bod.vJxS = vJxSb; bod.psidd = psiddb;
  bod.I = Ib;  bod.a = ab;% bod.BC = BCb;    
  bod.IC = ICb;  bod.f = fb;
  bod.Ia = Iab;
   
%   bod.X = 
%   bod.Sd = Sd;      bod.psid = psid;
%   bod.psidd = psidd;   bod.BC = BC;      bod.fC = f;
%   bod.v = v;           bod.a = a;        bod.vJ = vJ;  
%   bod.vJxS = vJxS;     bod.I = I;        bod.IC = IC;
%   bod.Ia = Ia;
%   
end



function Idot = dot(I, v)
    Idot = crf(v)*I - I*crm(v);
end

