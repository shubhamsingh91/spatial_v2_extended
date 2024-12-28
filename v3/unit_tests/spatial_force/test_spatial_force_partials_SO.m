clc; clear all;

% run([pwd,'\..\startup.m'])
% Testing the SO partials of spatial force

N = 17;

% Create a random model with N links
model = autoTree(N, 1.5, pi/3);
model.jtype{1}='Fb';
model = postProcessModel(model);

q  = rand(model.NQ,1);

if strcmp(model.jtype{1},'SE3')
    q(1:16) = reshape( randomSE3(),[16 1]);
elseif strcmp(model.jtype{1},'Fb')  % normalizing the quoternion
    q(1:4) = q(1:4)/norm(q(1:4));
end

qd = rand(model.NV,1);
qdd = rand(model.NV,1);

fprintf("N = %d\n", N)% 

fprintf("model.NQ = %d\n", model.NQ)% 
fprintf("model.NV = %d\n", model.NV)% 

%% complex-step- SO derivatives of spatial force/spatial cumulative force
step=1e-20; % step size for complex-step

newConfig = @(x) configurationAddition(model,q,x);
[glob,bod] = GlobalDynamics_red( model, q, qd, qdd);

for ii=1:N
    for jj=1:N
         for kk=1:N
                                                                                                                               % i variables
        Si = glob.S{ii};       psidi = glob.psid{ii};      fCi = glob.fC{ii}; 
        Sdi = glob.Sd{ii};      psiddi = glob.psidd{ii};
        ICi = glob.IC{ii};      BCi = glob.BC{ii};   
        BIic_Si = Bten(ICi,Si);
        BIic_psidi = Bten(ICi,psidi);
       
                                                                                                                  % j variables
        Sj = glob.S{jj};       psidj = glob.psid{jj};  psiddj = glob.psidd{jj};  
        Sdj = glob.Sd{jj};      
        
                                                                                                                  % k variables
        Sk = glob.S{kk};       psidk = glob.psid{kk}; 
        Sdk = glob.Sd{kk};     psiddk = glob.psidd{kk};
        
        BIic_psidj = Bten(ICi,psidj);
        
            if (ismember(kk,model.ancestors{jj})&&ismember(jj,model.ancestors{ii}))                                   % k<=j<=i
                
               fprintf("\n ii = % d; jj= %d; kk= %d  \n \n",ii,jj,kk)      
                   
               [d2fi_dqj_dqk_cs] = complexStepForce(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, ii , jj), ...
                zeros(model.NV,1), jj, kk);
                                         
                temp1 = 2*Tm(BIic_psidj + Tm(cmfM(Sj),BCi) - mT(BCi,crmM(Sj)) ,psidk) + ...
                            Tm( Tm(cmfM(Sj),ICi)- mT(ICi,crmM(Sj)) ,psiddk);
                
                temp2 = Tm(cmfM(Sk), 2*BCi*psidj + ICi*psiddj+ cmf_bar(fCi)*Sj);
                 
                d2fi_dqj_dqk = rotR(temp1) + temp2;
                
                if kk~=jj                                                                                             % k<j<= i
          
                   %------------------------ 
                   
                 [d2fi_dqk_dqj_cs] = complexStepForce(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, ii , kk), ...
                  zeros(model.NV,1), kk, jj);
                  d2fi_dqk_dqj = rotR(d2fi_dqj_dqk);
                   compare('(d2fic_dqk_dqj)'  , d2fi_dqk_dqj , d2fi_dqk_dqj_cs);
          
                   %------------------------
                 [d2fk_dqi_dqj_cs] = complexStepForce(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, kk , ii), ...
                  zeros(model.NV,1), ii, jj);
                   
                    temp1 = 2*Tm(BIic_psidi + Tm(cmfM(Si),BCi) - mT(BCi,crmM(Si)) ,psidj) + ...
                                Tm( Tm(cmfM(Si),ICi)- mT(ICi,crmM(Si)) ,psiddj);

                    temp2 = Tm(cmfM(Sj), 2*BCi*psidi + ICi*psiddi+ cmf_bar(fCi)*Si);

                    d2fk_dqi_dqj = rotR(temp1) + temp2;
                    compare('(d2fk_dqi_dqj)'  , d2fk_dqi_dqj , d2fk_dqi_dqj_cs);

                    %------------------------
                    if jj~=ii   % k < j < i 
                        
                    [d2fk_dqj_dqi_cs] = complexStepForce(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, kk , jj), ...
                        zeros(model.NV,1), jj, ii);                    
                      
                    d2fk_dqj_dqi = rotR(d2fk_dqi_dqj);
        
                     compare('(d2fk_dqj_dqi)'  , d2fk_dqj_dqi , d2fk_dqj_dqi_cs);
             
                    
                    end
 
                end
                

            end
                        
            
         end
    end
end


    
%% Functions
function compare(txt, v1, v2)
    if (size(v1,3))>1&&(size(v2,3)>1) 
        e = tens_norm(v1-v2);
    else
        e = norm(v1-v2);
    end
    if e > 1e-7
        x = 'X';
        fprintf('%12s = %.3e  %s\n',txt,e,x); 
        error('ErrorID:SpecificError', 'Expression not correct');

    else
        fprintf('%12s = %.3e  \x2713\n%s\n',txt,e);         
    end
    
end




