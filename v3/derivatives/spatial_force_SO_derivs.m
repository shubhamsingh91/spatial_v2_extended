function  [derivs] = spatial_force_SO_derivs(model, q, qd, qdd)

% SO Partial Derivatives for Cumulative spatial-force for multi-DoF joints
% Contributors - Shubham Singh, singh281@utexas.edu 

crf_bar = @(x)[-skew(x(1:3)), -skew(x(4:6)) ; -skew(x(4:6)) ,zeros(3,3)]; % cross product matrix for force vectors with a bar over top

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
    Sd{i} = crm(v{i})*S{i};
    BC{i} = 2*factorFunctions(IC{i},v{i});
    f{i}  =  IC{i}*a{i} + crf(v{i})*IC{i}*v{i};
    d2fc_dq{i}  =  zeros(6,model.NV,model.NV);

end


for i = model.NB:-1:1
        
        BCi = BC{i};
        ICi = IC{i};
        fCi = f{i};
        
    for p=1:model.nv(i) % looping over each DoF of joint j

        S_p    = S{i}(:,p); Sd_p = Sd{i}(:,p); 
        psid_p = psid{i}(:,p); psidd_p = psidd{i}(:,p);
        
        Bic_phii     = 2*factorFunctions(IC{i} ,S_p   );
        Bic_psii_dot = 2*factorFunctions(IC{i} ,psid_p);
       
        
      ii = model.vinds{i};
      j = i;

      while j > 0
        jj = model.vinds{j}; 
        
          for t=1:model.nv(j)
              
            S_t    = S{j}(:,t);    Sd_t   = Sd{j}(:,t);
            psid_t = psid{j}(:,t); psidd_t = psidd{j}(:,t);
        
            Bic_psijt_dot = 2*factorFunctions(IC{i} ,psid_t);
    
               k = j;
                  while k > 0
                    kk = model.vinds{k};
                    
                     for r=1:model.nv(k) % kk <= j <= i
                         
                        S_r    = S{k}(:,r);    Sd_r   = Sd{k}(:,r);
                        psid_r = psid{k}(:,r); psidd_r = psidd{k}(:,r);
            
                          
                           d2fc_dq{ii(p)}(:,jj(t),kk(r)) = (Bic_psijt_dot +2*dot(BCi,S_t))*psid_r +...
                                                            dot(ICi,S_t)*psidd_r +...
                                                            crf(S_r)*(BCi*psid_t+ICi*psidd_t+crf_bar(fCi)*S_t);  


                            if (j~=i)  % kk <= j < i

                               
                                
                            end

                             if (k~=j) % kk < j <= i
                                 
                                 

                                  if(j~=i) % kk < j < i

                                                      

                                  else % kk < j = i
                                      
                                  end
                                  
                             else     % kk = j <= i             

                             end

                          end
                    k  = model.parent(k);
                  end

        end 
        j  = model.parent(j);
      end
        
        
        
    end  

    if model.parent(i) > 0
        p = model.parent(i);
        IC{p} = IC{p} + IC{i};
        BC{p} = BC{p} + BC{i};
        f{p}  = f{p}  + f{i};
    end
end
derivs.d2fc_dq=d2fc_dq;


end

function Idot = dot(I, v)
    Idot = crf(v)*I - I*crm(v);
end

