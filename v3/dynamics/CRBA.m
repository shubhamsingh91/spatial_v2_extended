function  [M] = CRBA(model, q )

if ~iscell(q) 
    [q] = confVecToCell(model,q);
end

    a_grav = get_gravity(model);
    IC = model.I;

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

    S{i} = inv(X_m_0{i})*S{i}; % S in ground frame

end

    M  =  zeros(model.NV,model.NV);


for i = model.NB:-1:1
    
  ii = model.vinds{i};       
        
  F{i} = IC{i}*S{i};
  
   j=i;
   
   while (j>0)
        jj = model.vinds{j};
        t1 = S{j}.'*F{i};
        M(jj,ii) = t1;
        M(ii,jj) =   t1.';
        
        j = model.parent(j);
   end
   
    
    if model.parent(i) > 0
        p = model.parent(i);
        IC{p} = IC{p} + IC{i};
    end
end


end


