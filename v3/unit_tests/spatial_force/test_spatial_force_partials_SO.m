clc; clear all;

% run([pwd,'\..\startup.m'])
% Testing the SO partials of spatial force

N = 12;

% Create a random model with N links
model = autoTree(N, 2, pi/3);

% model.jtype{1}='Fb';
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

%% Testing cumulative spatial force derivs using full single-DoF algo

[derivs] = spatial_force_SO_derivs(model, q, qd, qdd);

d2fc_dq_algo = derivs.d2fc_dq;
d2fc_dv_algo = derivs.d2fc_dv;
d2fc_daq_algo = derivs.d2fc_daq;
d2fc_dqa_algo = derivs.d2fc_dqa;
d2fc_dvq_algo = derivs.d2fc_dvq;
d2fc_dva_algo = derivs.d2fc_dva;
d2fc_dav_algo = derivs.d2fc_dav;

d2fc_dq_algo_b = derivs.d2fc_dq_b;
d2fc_dv_algo_b = derivs.d2fc_dv_b;
d2fc_dvq_algo_b = derivs.d2fc_dvq_b;
d2fc_daq_algo_b = derivs.d2fc_daq_b;

%% complex-step- SO derivatives of spatial force/spatial cumulative force
step=1e-20; % step size for complex-step

newConfig = @(x) configurationAddition(model,q,x);
[glob,bod] = GlobalDynamics_red( model, q, qd, qdd);
tic

for ii=1:N
    for jj=1:N
         for kk=1:N
             
             ii_vec = model.vinds{ii}; jj_vec = model.vinds{jj};
             kk_vec = model.vinds{kk};
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
   
        BIic_Sj = Bten(ICi,Sj);
        BIic_psidj = Bten(ICi,psidj);
        BIic_psidk = Bten(ICi,psidk);

            if (ismember(kk,model.ancestors{jj})&&ismember(jj,model.ancestors{ii}))                                   % k<=j<=i
                
               fprintf("\n ii = % d; jj= %d; kk= %d  \n \n",ii,jj,kk)      
                   
               [d2fi_dqj_dqk_cs,d2fi_dqj_dqk_cs_body] = complexStepForce(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, ii , jj), ...
                zeros(model.NV,1), jj, kk);
                                         
                temp1 = 2*Tm(BIic_psidj + Tm(cmfM(Sj),BCi) - mT(BCi,crmM(Sj)) ,psidk) + ...
                            Tm( Tm(cmfM(Sj),ICi)- mT(ICi,crmM(Sj)) ,psiddk);
                
                temp2 = Tm(cmfM(Sk), 2*BCi*psidj + ICi*psiddj+ cmf_bar(fCi)*Sj);
                 
                d2fi_dqj_dqk = rotR(temp1) + temp2;
                
                compare('(d2fic_dqj_dqk)'  , d2fi_dqj_dqk , d2fi_dqj_dqk_cs);
                compare('(d2fic_dqj_dqk-- algo)'  , d2fi_dqj_dqk , d2fc_dq_algo{ii}(1:6,jj_vec,kk_vec));
                compare('(d2fic_dqj_dqk) - body'  , d2fi_dqj_dqk_cs_body ,  d2fc_dq_algo_b{ii}(1:6,jj_vec,kk_vec));

                %------- SO a/q Case 1A
                
                [d2fi_daj_dqk_cs,d2fi_daj_dqk_cs_b] = complexStepForce_SOaq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, ii , jj), ...
                zeros(model.NV,1), jj, kk);
                                         
                d2fi_daj_dqk = Tm(Tm(cmfM(Sk),ICi) , Sj);
                
                compare('(d2fic_daj_dqk) case 1A'  , d2fi_daj_dqk , d2fi_daj_dqk_cs);               
                compare('(d2fic_daj_dqk) case 1A - algo'  , d2fi_daj_dqk , d2fc_daq_algo{ii}(:,jj_vec,kk_vec));               
                compare('(d2fic_daj_dqk) case 1A - algo --body'  , d2fi_daj_dqk_cs_b , d2fc_daq_algo_b{ii}(:,jj_vec,kk_vec));               
              
                %------- SO q/a Case 1A
                [d2fi_dqj_dak_cs] = complexStepForce_SOqa(model, @(x) spatial_force_derivatives(model, q ,qd ,x, ii , jj), ...
                zeros(model.NV,1), jj, kk);
                compare('(d2fi_dqj_dak_cs) case 1A - algo'  , d2fi_dqj_dak_cs , d2fc_dqa_algo{ii}(:,jj_vec,kk_vec));               

                %------- SO v/q Case 1A
                
                [d2fi_dvj_dqk_cs,d2fi_dvj_dqk_cs_b] = complexStepForce_SOvq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, ii , jj), ...
                zeros(model.NV,1), jj, kk);
                                         
                 
                d2fi_dvj_dqk = 2*Tm(BIic_psidk + Tm(cmfM(Sk),BCi)+mT(ICi,crmM(psidk)), Sj) +...
                            Tm(cmfM(Sk),ICi*(psidj+Sdj));
                
                compare('(d2fi_dvj_dqk) case 1A'  , d2fi_dvj_dqk , d2fi_dvj_dqk_cs);
                compare('(d2fi_dvj_dqk) case 1A -- algo'  , d2fi_dvj_dqk , d2fc_dvq_algo{ii}(:,jj_vec,kk_vec));
                compare('(d2fi_dvj_dqk) case 1A -- algo --body'  , d2fi_dvj_dqk_cs_b , d2fc_dvq_algo_b{ii}(:,jj_vec,kk_vec));

                if kk~=jj  % k<j<=i                                                                                           % k<j<= i
          
                   %------------------------ 
                   
                 [d2fi_dqk_dqj_cs, d2fi_dqk_dqj_cs_bod] = complexStepForce(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, ii , kk), ...
                  zeros(model.NV,1), kk, jj);
                  d2fi_dqk_dqj = rotR(d2fi_dqj_dqk);
                   compare('(d2fic_dqk_dqj)'  , d2fi_dqk_dqj , d2fi_dqk_dqj_cs);
                   compare('(d2fic_dqk_dqj) - algo'  , d2fi_dqk_dqj ,  d2fc_dq_algo{ii}(1:6,kk_vec,jj_vec));
                   compare('(d2fic_dqk_dqj) - algo -- body'  , d2fi_dqk_dqj_cs_bod ,  d2fc_dq_algo_b{ii}(1:6,kk_vec,jj_vec));

                   %------------------------
                 [d2fk_dqi_dqj_cs,d2fk_dqi_dqj_cs_bod] = complexStepForce(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, kk , ii), ...
                  zeros(model.NV,1), ii, jj);
                   
                    temp1 = 2*Tm(BIic_psidi + Tm(cmfM(Si),BCi) - mT(BCi,crmM(Si)) ,psidj) + ...
                                Tm( Tm(cmfM(Si),ICi)- mT(ICi,crmM(Si)) ,psiddj);

                    temp2 = Tm(cmfM(Sj), 2*BCi*psidi + ICi*psiddi+ cmf_bar(fCi)*Si);

                    d2fk_dqi_dqj = rotR(temp1) + temp2;
                    compare('(d2fk_dqi_dqj)'  , d2fk_dqi_dqj , d2fk_dqi_dqj_cs);
                    compare('(d2fk_dqi_dqj- algo)'  , d2fk_dqi_dqj ,  d2fc_dq_algo{kk}(1:6,ii_vec,jj_vec));
                    compare('(d2fk_dqi_dqj- algo) --body'  , d2fk_dqi_dqj_cs_bod ,  d2fc_dq_algo_b{kk}(1:6,ii_vec,jj_vec));

                     
                     %------- SO a/q Case 2B j !=i

                    [d2fk_dai_dqj_cs,d2fk_dai_dqj_cs_b] = complexStepForce_SOaq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, kk , ii), ...
                    zeros(model.NV,1), ii, jj);


                    d2fk_dai_dqj = Tm(cmfM(Sj), ICi*Si);

                    compare('(d2fk_dai_dqj) case 2B'  , d2fk_dai_dqj , d2fk_dai_dqj_cs); 
                    compare('(d2fk_dai_dqj) case 2B --algo'  , d2fk_dai_dqj , d2fc_daq_algo{kk}(:,ii_vec,jj_vec));                    
                    compare('(d2fk_dai_dqj) case 2B --algo --body'  , d2fk_dai_dqj_cs_b , d2fc_daq_algo_b{kk}(:,ii_vec,jj_vec));                    

          
                    %--te----- SO v/q Case 2B

                    [d2fk_dvi_dqj_cs,d2fk_dvi_dqj_cs_b] = complexStepForce_SOvq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, kk , ii), ...
                    zeros(model.NV,1), ii, jj);

                    d2fk_dvi_dqj = rotR(Tm(2*BIic_Si,psidj)+ Tm(cmf_barM(2*BCi*Si + ICi*(psidi+Sdi)),Sj));

                    compare('(d2fk_dvi_dqj) case 2B'  , d2fk_dvi_dqj , d2fk_dvi_dqj_cs);
                    compare('(d2fk_dvi_dqj) case 2B-- algo'  , d2fk_dvi_dqj , d2fc_dvq_algo{kk}(:,ii_vec,jj_vec));
                    compare('(d2fk_dvi_dqj) case 2B-- algo--body'  , d2fk_dvi_dqj_cs_b , d2fc_dvq_algo_b{kk}(:,ii_vec,jj_vec));

                    %------------------------
                    if jj~=ii   % k < j < i 
                        
                    [d2fk_dqj_dqi_cs] = complexStepForce(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, kk , jj), ...
                        zeros(model.NV,1), jj, ii);                    
                      
                    d2fk_dqj_dqi = rotR(d2fk_dqi_dqj);
        
                     compare('(d2fk_dqj_dqi)'  , d2fk_dqj_dqi , d2fk_dqj_dqi_cs);
                     compare('(d2fk_dqj_dqi) - algo'  , d2fk_dqj_dqi , d2fc_dq_algo{kk}(1:6,jj_vec,ii_vec));
                     
                     
                     %------- SO w.r.t v for Case 1C
                    [d2fk_dvi_dvj_cs,d2fk_dvi_dvj_cs_bod] = complexStepForce_SOv(model, @(x) spatial_force_derivatives(model, q ,x ,qdd, kk , ii), ...
                    zeros(model.NV,1), ii, jj);
                
                   d2fk_dvi_dvj = rotR(Tm(2*BIic_Si,Sj));
               
                  compare('(d2fk_dvi_dvj) case 1C'  , d2fk_dvi_dvj , d2fk_dvi_dvj_cs);                
                  compare('(d2fk_dvi_dvj) case 1C-- algo'  , d2fk_dvi_dvj , d2fc_dv_algo{kk}(:,ii_vec,jj_vec));                
                  compare('(d2fk_dvi_dvj) case 1C-- algo- body'  , d2fk_dvi_dvj_cs_bod , d2fc_dv_algo_b{kk}(:,ii_vec,jj_vec));                

                       %------- SO w.r.t v for Case 2C
                     [d2fk_dvj_dvi_cs,d2fk_dvj_dvi_cs_body] = complexStepForce_SOv(model, @(x) spatial_force_derivatives(model, q ,x ,qdd, kk , jj), ...
                    zeros(model.NV,1), jj, ii);
                
                    d2fk_dvj_dvi = rotR(d2fk_dvi_dvj);
               
                    compare('(d2fk_dvj_dvi) case 2C'  , d2fk_dvj_dvi , d2fk_dvj_dvi_cs);           
                    compare('(d2fk_dvj_dvi) case 2C -- algo'  , d2fk_dvj_dvi , d2fc_dv_algo{kk}(:,jj_vec,ii_vec));           
                    compare('(d2fk_dvj_dvi) case 2C -- algo --body'  , d2fk_dvj_dvi_cs_body , d2fc_dv_algo_b{kk}(:,jj_vec,ii_vec));           

                %------- SO q/a Case 2B j !=i
                [d2fk_dqi_daj_cs] = complexStepForce_SOqa(model, @(x) spatial_force_derivatives(model, q ,qd ,x, kk , ii), ...
                zeros(model.NV,1), ii, jj);
               
                compare('(d2fk_dqi_daj) case 2B --algo'  , d2fk_dqi_daj_cs , d2fc_dqa_algo{kk}(:,ii_vec,jj_vec)); 

               %------- SO a/q Case 2C
                
                [d2fk_daj_dqi_cs,d2fk_daj_dqi_cs_body] = complexStepForce_SOaq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, kk , jj), ...
                zeros(model.NV,1), jj, ii);
                                         
                d2fk_daj_dqi = Tm(Tm(cmfM(Si),ICi) - mT(ICi,crmM(Si)), Sj);
                
                compare('(d2fk_daj_dqi) Case 2C'  , d2fk_daj_dqi , d2fk_daj_dqi_cs); 
                compare('(d2fk_daj_dqi) Case 2C -- algo'  , d2fk_daj_dqi , d2fc_daq_algo{kk}(:,jj_vec,ii_vec)); 
                compare('(d2fk_daj_dqi) Case 2C -- algo --body'  , d2fk_daj_dqi_cs_body , d2fc_daq_algo_b{kk}(:,jj_vec,ii_vec)); 

                %------- SO q/a Case 2C
                [d2fk_dqj_dai_cs] = complexStepForce_SOqa(model, @(x) spatial_force_derivatives(model, q ,qd ,x, kk , jj), ...
                zeros(model.NV,1), jj, ii);
            
                compare('(d2fk_dqj_dai) Case 2C -- algo'  , d2fk_dqj_dai_cs , d2fc_dqa_algo{kk}(:,jj_vec,ii_vec)); 


                  %------- SO v/q Case 2C
                
                [d2fk_dvj_dqi_cs,d2fk_dvj_dqi_cs_b] = complexStepForce_SOvq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, kk , jj), ...
                zeros(model.NV,1), jj, ii);
                                         
                d2fk_dvj_dqi = 2*Tm(BIic_psidi + Tm(cmfM(Si),BCi)- mT(BCi,crmM(Si)), Sj) +...
                            Tm(Tm(cmfM(Si),ICi) - mT(ICi,crmM(Si))  ,psidj+Sdj);
                
                compare('(d2fk_dvj_dqi) case 2C'  , d2fk_dvj_dqi , d2fk_dvj_dqi_cs);
                compare('(d2fk_dvj_dqi) case 2C -- algo'  , d2fk_dvj_dqi , d2fc_dvq_algo{kk}(:,jj_vec,ii_vec));
                compare('(d2fk_dvj_dqi) case 2C -- algo --body'  , d2fk_dvj_dqi_cs_b , d2fc_dvq_algo_b{kk}(:,jj_vec,ii_vec));

                    else % k < j = i
                        
                     [d2fk_dvi_dvj_cs,d2fk_dvi_dvj_cs_body] = complexStepForce_SOv(model, @(x) spatial_force_derivatives(model, q ,x ,qdd, kk , ii), ...
                    zeros(model.NV,1), ii, jj);
                
                    d2fk_dvi_dvj = rotR(Tm (cmf_barM(ICi*Si)+ Tm(cmfM(Si),ICi) ,Sj));
               
                    compare('(d2fk_dvi_dvj) case D'  , d2fk_dvi_dvj , d2fk_dvi_dvj_cs);     
                    compare('(d2fk_dvi_dvj) case D -- algo'  , d2fk_dvi_dvj , d2fc_dv_algo{kk}(:,ii_vec,jj_vec));     
                    compare('(d2fk_dvi_dvj) case D -- algo --body'  , d2fk_dvi_dvj_cs_body , d2fc_dv_algo_b{kk}(:,ii_vec,jj_vec));     
                       
                    end
 
                %--------- SO derivs w.r.t v case 1A
                      
               [d2fi_dvj_dvk_cs,d2fi_dvj_dvk_cs_body] = complexStepForce_SOv(model, @(x) spatial_force_derivatives(model, q ,x ,qdd, ii , jj), ...
                zeros(model.NV,1), jj, kk);
                
                d2fi_dvj_dvk = rotR(Tm(2*BIic_Sj,Sk));
               
                compare('(d2fi_dvj_dvk)'  , d2fi_dvj_dvk , d2fi_dvj_dvk_cs);
                compare('(d2fi_dvj_dvk) - algo'  , d2fi_dvj_dvk , d2fc_dv_algo{ii}(:,jj_vec,kk_vec));
                compare('(d2fi_dvj_dvk) - algo --body'  , d2fi_dvj_dvk_cs_body , d2fc_dv_algo_b{ii}(:,jj_vec,kk_vec));
                
                %--------- SO derivs w.r.t v case 2A
                      
               [d2fi_dvk_dvj_cs,d2fi_dvk_dvj_cs_body] = complexStepForce_SOv(model, @(x) spatial_force_derivatives(model, q ,x ,qdd, ii , kk), ...
                zeros(model.NV,1), kk, jj);
                
                d2fi_dvk_dvj = rotR(d2fi_dvj_dvk);
               
                compare('(d2fi_dvk_dvj)'  , d2fi_dvk_dvj , d2fi_dvk_dvj_cs);
                compare('(d2fi_dvk_dvj) -- algo'  , d2fi_dvk_dvj , d2fc_dv_algo{ii}(:,kk_vec,jj_vec));
                compare('(d2fi_dvk_dvj) -- algo -body'  , d2fi_dvk_dvj_cs_body , d2fc_dv_algo_b{ii}(:,kk_vec,jj_vec));

                %------- SO a/q Case 1B
                
                [d2fi_dak_dqj_cs,d2fi_dak_dqj_cs_b] = complexStepForce_SOaq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, ii , kk), ...
                zeros(model.NV,1), kk, jj);
                                         
                d2fi_dak_dqj = Tm(Tm(cmfM(Sj),ICi) - mT(ICi,crmM(Sj)) , Sk);
                
                compare('(d2fi_dak_dqj) case 1A'  , d2fi_dak_dqj , d2fi_dak_dqj_cs);
                compare('(d2fi_dak_dqj) case 1A -- algo'  , d2fi_dak_dqj , d2fc_daq_algo{ii}(:,kk_vec,jj_vec));
                compare('(d2fi_dak_dqj) case 1A -- algo --body'  , d2fi_dak_dqj_cs_b , d2fc_daq_algo_b{ii}(:,kk_vec,jj_vec));

               %------- SO q/a Case 1B
                [d2fi_dqk_daj_cs] = complexStepForce_SOqa(model, @(x) spatial_force_derivatives(model, q ,qd ,x, ii , kk), ...
                zeros(model.NV,1), kk, jj);
            
                compare('(d2fi_dqk_daj) case 1A -- algo'  , d2fi_dqk_daj_cs , d2fc_dqa_algo{ii}(:,kk_vec,jj_vec));

               %------- SO v/q Case 1B
                
                [d2fi_dvk_dqj_cs,d2fi_dvk_dqj_cs_b] = complexStepForce_SOvq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, ii , kk), ...
                zeros(model.NV,1), kk, jj);
                                         
                d2fi_dvk_dqj = 2*Tm(BIic_psidj + Tm(cmfM(Sj),BCi)- mT(BCi,crmM(Sj)), Sk) +...
                            Tm(Tm(cmfM(Sj),ICi) - mT(ICi,crmM(Sj))  ,psidk+Sdk);
                
                compare('(d2fi_dvk_dqj) case 1B'  , d2fi_dvk_dqj , d2fi_dvk_dqj_cs);
                compare('(d2fi_dvk_dqj) case 1B -- algo'  , d2fi_dvk_dqj , d2fc_dvq_algo{ii}(:,kk_vec,jj_vec));
                compare('(d2fi_dvk_dqj) case 1B -- algo --body'  , d2fi_dvk_dqj_cs_b , d2fc_dvq_algo_b{ii}(:,kk_vec,jj_vec));

                
                else % k=j<=i
                 %--------- SO derivs w.r.t v case B
                 
               [d2fi_dvj_dvk_cs,d2fi_dvj_dvk_cs_body] = complexStepForce_SOv(model, @(x) spatial_force_derivatives(model, q ,x ,qdd, ii , jj), ...
                zeros(model.NV,1), jj, kk);
            
                d2fi_dvj_dvk = rotR(Tm( Tm(cmfM(Sj),ICi)+cmf_barM(ICi*Sj), Sk));
                  compare('(d2fi_dvk_dvj) Case B'  , d2fi_dvj_dvk , d2fi_dvj_dvk_cs);
                  compare('(d2fi_dvk_dvj) Case B -- algo'  , d2fi_dvj_dvk , d2fc_dv_algo{ii}(:,jj_vec,kk_vec));
                  compare('(d2fi_dvk_dvj) Case B -- algo -- body'  , d2fi_dvj_dvk_cs_body , d2fc_dv_algo_b{ii}(:,jj_vec,kk_vec));
                 
                end
                
                %-----
                if jj~= ii % kk <= jj < ii
                    
                [d2fj_dqk_dqi_cs,d2fj_dqk_dqi_cs_body] = complexStepForce(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, jj , kk), ...
                        zeros(model.NV,1), kk, ii);
                    
                    
                 temp1 = 2*Tm(BIic_psidi + Tm(cmfM(Si),BCi) - mT(BCi,crmM(Si)) ,psidk) + ...
                                Tm( Tm(cmfM(Si),ICi)- mT(ICi,crmM(Si)) ,psiddk);

                 temp2 = Tm(cmf_barM(2*BCi*psidi + ICi*psiddi+ cmf_bar(fCi)*Si), Sk);

                 d2fj_dqk_dqi = temp1 + temp2;                   
                      compare('(d2fj_dqk_dqi)'  , d2fj_dqk_dqi , d2fj_dqk_dqi_cs);
                      compare('(d2fj_dqk_dqi)-- algo'  , d2fj_dqk_dqi , d2fc_dq_algo{jj}(1:6,kk_vec,ii_vec));
                      compare('(d2fj_dqk_dqi)-- algo --body'  , d2fj_dqk_dqi_cs_body , d2fc_dq_algo_b{jj}(1:6,kk_vec,ii_vec));

                %-----
                [d2fj_dqi_dqk_cs,d2fj_dqi_dqk_cs_bod] = complexStepForce(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, jj , ii), ...
                        zeros(model.NV,1), ii, kk);
                    
                d2fj_dqi_dqk = rotR(d2fj_dqk_dqi);
                
                 compare('(d2fj_dqi_dqk)'  , d2fj_dqi_dqk , d2fj_dqi_dqk_cs);
                 compare('(d2fj_dqi_dqk) -- algo'  , d2fj_dqi_dqk , d2fc_dq_algo{jj}(1:6,ii_vec,kk_vec));
                 compare('(d2fj_dqi_dqk) -- algo-- body'  , d2fj_dqi_dqk_cs_bod , d2fc_dq_algo_b{jj}(1:6,ii_vec,kk_vec));

                      
                %----- SO v case 1E
                [d2fj_dvk_dvi_cs, d2fj_dvk_dvi_cs_body] = complexStepForce_SOv(model, @(x) spatial_force_derivatives(model, q ,x ,qdd, jj , kk), ...
                zeros(model.NV,1), kk, ii);
                
                d2fj_dvk_dvi = 2*Tm(BIic_Si,Sk);
               
                compare('(d2fj_dvk_dvi)'  , d2fj_dvk_dvi , d2fj_dvk_dvi_cs);
                compare('(d2fj_dvk_dvi) -- algo'  , d2fj_dvk_dvi , d2fc_dv_algo{jj}(:,kk_vec,ii_vec));
                compare('(d2fj_dvk_dvi) -- algo -- body'  , d2fj_dvk_dvi_cs_body , d2fc_dv_algo_b{jj}(:,kk_vec,ii_vec));

                %----- SO v case 2E
                [d2fj_dvi_dvk_cs, d2fj_dvi_dvk_cs_body] = complexStepForce_SOv(model, @(x) spatial_force_derivatives(model, q ,x ,qdd, jj , ii), ...
                zeros(model.NV,1), ii, kk);
                
                d2fj_dvi_dvk = rotR(d2fj_dvk_dvi);
               
                compare('(d2fj_dvi_dvk)'  , d2fj_dvi_dvk , d2fj_dvi_dvk_cs);
                compare('(d2fj_dvi_dvk) -- algo'  , d2fj_dvi_dvk , d2fc_dv_algo{jj}(:,ii_vec,kk_vec));
                compare('(d2fj_dvi_dvk) -- algo -- body'  , d2fj_dvi_dvk_cs_body , d2fc_dv_algo_b{jj}(:,ii_vec,kk_vec));
                
                %------- SO a/q Case 2A
                
                [d2fj_dai_dqk_cs,d2fj_dai_dqk_cs_b] = complexStepForce_SOaq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, jj , ii), ...
                zeros(model.NV,1), ii, kk);
                                         
                d2fj_dai_dqk = Tm(Tm(cmfM(Sk),ICi) , Si);
                
                compare('(d2fj_dai_dqk) Case 2A'  , d2fj_dai_dqk , d2fj_dai_dqk_cs);               
                compare('(d2fj_dai_dqk) Case 2A -- algo'  , d2fj_dai_dqk , d2fc_daq_algo{jj}(:,ii_vec,kk_vec));               
                compare('(d2fj_dai_dqk) Case 2A -- algo --body'  , d2fj_dai_dqk_cs_b , d2fc_daq_algo_b{jj}(:,ii_vec,kk_vec));               
             
                %------- SO q/a Case 2A
                [d2fj_dqi_dak_cs] = complexStepForce_SOqa(model, @(x) spatial_force_derivatives(model, q ,qd ,x, jj , ii), ...
                zeros(model.NV,1), ii, kk);

                compare('(d2fj_dqi_dak) Case 2A -- algo'  , d2fj_dqi_dak_cs , d2fc_dqa_algo{jj}(:,ii_vec,kk_vec));               
           
                 
               %------- SO a/q Case 1C
                
                [d2fj_dak_dqi_cs,d2fj_dak_dqi_cs_b] = complexStepForce_SOaq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, jj , kk), ...
                zeros(model.NV,1), kk, ii);
                                         
                d2fj_dak_dqi = Tm(Tm(cmfM(Si),ICi) - mT(ICi,crmM(Si)), Sk);
                
                compare('(d2fj_dak_dqi) Case 2A'  , d2fj_dak_dqi , d2fj_dak_dqi_cs); 
                compare('(d2fj_dak_dqi) Case 2A-- algo'  , d2fj_dak_dqi , d2fc_daq_algo{jj}(:,kk_vec,ii_vec)); 
                compare('(d2fj_dak_dqi) Case 2A-- algo --body'  , d2fj_dak_dqi_cs_b , d2fc_daq_algo_b{jj}(:,kk_vec,ii_vec)); 
              
               %------- SO q/a Case 1C
                [d2fj_dqk_dai_cs] = complexStepForce_SOqa(model, @(x) spatial_force_derivatives(model,q ,qd ,x, jj , kk), ...
                zeros(model.NV,1), kk, ii);
               
              compare('(d2fj_dqk_dai) Case 2A-- algo'  , d2fj_dqk_dai_cs , d2fc_dqa_algo{jj}(:,kk_vec,ii_vec)); 

            
              %------- SO v/q Case 2A
                
                [d2fj_jdvi_dqk_cs,d2fj_jdvi_dqk_cs_b] = complexStepForce_SOvq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, jj , ii), ...
                zeros(model.NV,1), ii, kk);
                                         
                 
                d2fj_dvi_dqk = 2*Tm(BIic_psidk + Tm(cmfM(Sk),BCi)+mT(ICi,crmM(psidk)), Si) +...
                            Tm(cmfM(Sk),ICi*(psidi+Sdi));
                
                compare('(d2fj_dvi_dqk) case 2A'  , d2fj_dvi_dqk , d2fj_jdvi_dqk_cs);
                compare('(d2fj_dvi_dqk) case 2A -- algo'  , d2fj_dvi_dqk , d2fc_dvq_algo{jj}(:,ii_vec,kk_vec));
                compare('(d2fj_dvi_dqk) case 2A -- algo --body'  , d2fj_jdvi_dqk_cs_b , d2fc_dvq_algo_b{jj}(:,ii_vec,kk_vec));

               %------- SO v/q Case 1C
                
                [d2fj_dvk_dqi_cs,d2fj_dvk_dqi_cs_b] = complexStepForce_SOvq(model, @(x) spatial_force_derivatives(model, newConfig(x) ,qd ,qdd, jj , kk), ...
                zeros(model.NV,1), kk, ii);
                                         
                d2fj_dvk_dqi = 2*Tm(BIic_psidi + Tm(cmfM(Si),BCi)- mT(BCi,crmM(Si)), Sk) +...
                            Tm(Tm(cmfM(Si),ICi) - mT(ICi,crmM(Si))  ,psidk+Sdk);
                
                compare('(d2fj_dvk_dqi) case 1C'  , d2fj_dvk_dqi , d2fj_dvk_dqi_cs);
                compare('(d2fj_dvk_dqi) case 1C -- algo'  , d2fj_dvk_dqi , d2fc_dvq_algo{jj}(:,kk_vec,ii_vec));
                compare('(d2fj_dvk_dqi) case 1C -- algo --body'  , d2fj_dvk_dqi_cs_b , d2fc_dvq_algo_b{jj}(:,kk_vec,ii_vec));

                end
                
            % SO-a 
               [d2fi_daj_dak_cs] = complexStepForce_SOa(model, @(x) spatial_force_derivatives(model, q ,qd ,x, ii , jj), ...
                zeros(model.NV,1), jj, kk);
                
                d2fi_daj_dak = zeros(size(d2fi_daj_dak_cs));
                 compare('(d2fi_daj_dak) -- all cases'  , d2fi_daj_dak , d2fi_daj_dak_cs);
               
                
            % SO-a/v 
               [d2fi_daj_dvk_cs] = complexStepForce_SOav(model, @(x) spatial_force_derivatives(model, q ,x ,qdd, ii , jj), ...
                zeros(model.NV,1), jj, kk);
                
                d2fi_daj_dvk = zeros(size(d2fi_daj_dvk_cs));
                 compare('(d2fi_daj_dvk) -- all cases'  , d2fi_daj_dvk , d2fi_daj_dvk_cs);
                 compare('(d2fi_daj_dvk) -- algo'  , d2fi_daj_dvk , d2fc_dav_algo{ii}(:,jj_vec,kk_vec));

            % SO-v/a
               [d2fi_dvj_dak_cs] = complexStepForce_SOva(model, @(x) spatial_force_derivatives(model, q ,qd ,x, ii , jj), ...
                zeros(model.NV,1), jj, kk);
                
                d2fi_dvj_dak = zeros(size(d2fi_dvj_dak_cs));
                 compare('(d2fi_dvj_dak) -- all cases'  , d2fi_dvj_dak , d2fi_dvj_dak_cs);               
                 compare('(d2fi_dvj_dak) -- algo'  , d2fi_dvj_dak , d2fc_dva_algo{ii}(:,jj_vec,kk_vec));               

            end

                        
            
         end
    end
end


toc
    
%% Functions
function compare(txt, v1, v2)

    if (size(v1)~= size(v2))
        error('ErrorID:SpecificError', 'Size not correct');
    end
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




