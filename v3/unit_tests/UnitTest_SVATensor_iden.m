% Unit tests for SVA Tensor identites 

clear
N = 5;

% Create a random model with N links
model = autoTree(N, 1.5, pi/3);
checkSVATensor_iden(model,'Floating Base No Rotors');

function checkSVATensor_iden(model, desc)
fprintf('====================================\n');
fprintf('%s\n',desc);
fprintf('====================================\n');

model.jtype{1}='Fb';
model = postProcessModel(model);

q  = rand(model.NQ,1);

if strcmp(model.jtype{1},'SE3')
    q(1:16) = reshape( randomSE3(),[16 1]);
elseif strcmp(model.jtype{1},'Fb')  % normalizing the quoternion
    q(1:4) = q(1:4)/norm(q(1:4));
end

qd = rand(model.NV,1);
qdd=rand(model.NV,1);
m6 = rand(6,1); % random 6-motion vector
f6 = rand(6,1); %random 6-force vector

%%
step=1e-20; % step size for complex-step

newConfig = @(x) configurationAddition(model,q,x);
[glob,bod] = GlobalDynamics( model, q, qd, qdd,m6,f6);
% For all possible ii and jj cases

for ii=1:model.NB                           
     for jj=1:ii

        ni = model.nv(ii); nj = model.nv(jj);
        p = randi([1,nj]); % 1 random single mode of joint-j

        [glob_cs_q] = perturbIndex_q(model,@(x) GlobalDynamics(model, newConfig(x) ...
            ,qd ,qdd,m6,f6), zeros(model.NV,1),jj,ii,p );                                                        % complex with perturb q

        glob_cs_v = perturbIndex_v(model,@(x) GlobalDynamics(model, q ...
            ,x ,qdd,m6,f6), qd,jj,ii );                                                                        % complex with perturb v

                                                                                                         % glob coordinates
        Si_g = glob.S{ii};       psidi_g = glob.psid{ii};  vi_g = glob.v{ii};      ai_g = glob.a{ii};
        Sdi_g = glob.Sd{ii};     Ii_g = glob.I{ii};        Iai_g = glob.Ia{ii};    psiddi_g = glob.psidd{ii};
        vJxSi_g = glob.vJxS{ii}; ICi_g = glob.IC{ii};      BCi_g = glob.BC{ii};    fCi_g = glob.fC{ii}; 
        fi_g = glob.f{ii};        Bi_g = glob.B{ii}; 
        
        Sj_g = glob.S{jj};       psidj_g = glob.psid{jj};  psiddj_g = glob.psidd{jj};  
        Sdj_g = glob.Sd{jj};     ICj_g = glob.IC{jj};      BCj_g = glob.BC{jj};
        fCj_g = glob.fC{jj};
        BIic_psidj_g = Bten(ICi_g,psidj_g);                BIjc_psidj_g = Bten(ICj_g,psidj_g);
        BIic_Sj_g = Bten(ICi_g,Sj_g);                      BIjc_Sj_g = Bten(ICj_g,Sj_g);
           
        
        dSi_dqj_g = zeros(6,ni,nj);             
        dSdi_dqj_g = zeros(6,ni,nj);             
        dvJxsi_dqj_g = zeros(6,ni,nj);         
        dpsidi_dqj_g = zeros(6,ni,nj);           
        dIi_dqj_g = zeros(6,6,nj);               
        dICi_dqj_g = zeros(6,6,nj);              
        dai_dqj_g = zeros(6,nj);                
        dIai_dqj_g = zeros(6,nj);     
        dpsiddi_dqj_g = zeros(6,ni,nj);   
        dBCi_dqj_g = zeros(6,6,nj);  
        dfi_dqj_g = zeros(6,nj);           
        dfCi_dqj_g = zeros(6,nj);      
        dSiT_dqj_g = zeros(ni,6,nj);     
        dvi_dqj_g = zeros(6,nj);
        
        dSdi_dvj_g = zeros(6,ni,nj);  
        dpsidi_dvj_g = zeros(6,ni,nj);  
        dBCi_dvj_g = zeros(6,6,nj);  
    
        if  ismember(jj,model.ancestors{ii})                                                             % jj<=ii

            dSi_dqj_g = Tm(crmM(Sj_g),Si_g);                                                             % K1
            dSdi_dqj_g = Tm(crmM(psidj_g),Si_g)+Tm(crmM(Sj_g),Sdi_g);                                    % K2
            dvJxsi_dqj_g = Tm(crmM(Sj_g),vJxSi_g);                                                       % K3
            dpsidi_dqj_g = Tm(crmM(psidj_g),Si_g)+Tm(crmM(Sj_g),psidi_g);                                % K4
            dIi_dqj_g = Tm(cmfM(Sj_g),Ii_g)-mT(Ii_g,crmM(Sj_g));                                         % K5
            dICi_dqj_g = Tm(cmfM(Sj_g),ICi_g)-mT(ICi_g,crmM(Sj_g));                                      % K6
            dai_dqj_g = psiddj_g - crm(vi_g)*psidj_g - crm(ai_g)*Sj_g;                                   % K7
            dIai_dqj_g = cmf_bar(Iai_g)*Sj_g + Ii_g*psiddj_g - Ii_g*(crm(vi_g)*psidj_g);                 % K8
            dpsiddi_dqj_g = Tm(crmM(psiddj_g),Si_g)+2*Tm(crmM(psidj_g),psidi_g)+Tm(crmM(Sj_g),psiddi_g); % K9
            dBCi_dqj_g = BIic_psidj_g+ Tm(cmfM(Sj_g),BCi_g) - mT(BCi_g,crmM(Sj_g));                      % K10
            dfi_dqj_g = Ii_g*psiddj_g+ cmf_bar(fi_g)*Sj_g + 2*Bi_g*psidj_g;                              % K11
            dfCi_dqj_g = ICi_g*psiddj_g + cmf_bar(fCi_g)*Sj_g + 2*BCi_g*psidj_g;                         % K12
            dSiT_dqj_g = -mT(Si_g.',cmfM(Sj_g));                                                         % K13
            dSdi_dvj_g = Tm(crmM(Sj_g),Si_g);                                                            % K14
            dBCi_dvj_g = BIic_Sj_g;                                                                      % K16
   

            if jj==ii
            fprintf("\n ii = %d; jj= %d  j<=i Case \n \n",ii,jj)
            else
              fprintf("\n ii = %d; jj= %d  j<i Case \n \n",ii,jj) 
               dpsidi_dvj_g = Tm(crmM(Sj_g),Si_g);                                                       % K15
            end


            elseif ismember(ii,model.ancestors{jj}) 
                 fprintf("\n ii = % d; jj= %d  j>i Case \n \n",ii,jj)  
                 dICi_dqj_g = Tm(cmfM(Sj_g),ICj_g)-mT(ICj_g,crmM(Sj_g));                                 % K6
                 dBCi_dqj_g = BIjc_psidj_g+ Tm(cmfM(Sj_g),BCj_g) - mT(BCj_g,crmM(Sj_g));                 % K10
                 dfCi_dqj_g = ICj_g*psiddj_g + cmf_bar(fCj_g)*Sj_g + 2*BCj_g*psidj_g;                    % K12
                 dBCi_dvj_g = BIjc_Sj_g;                                                                 % K16

            else
                  fprintf("\n ii = % d; jj= %d  different branches  \n \n",ii,jj)  

        end
        fprintf('Global Coordinates\n');

        compare('(K1) '  , dSi_dqj_g , glob_cs_q.dSi_dqj_cs);
        compare('(K2) '  , dSdi_dqj_g , glob_cs_q.dSdi_dqj_cs);
        compare('(K3) '  , dvJxsi_dqj_g ,  glob_cs_q.dvJxsi_dqj_cs);
        compare('(K4) '  , dpsidi_dqj_g ,  glob_cs_q.dpsidi_dqj_cs);
        compare('(K5) '  , dIi_dqj_g ,  glob_cs_q.dIi_dqj_cs);
        compare('(K6) '  , dICi_dqj_g ,  glob_cs_q.dICi_dqj_cs);
        compare('(K7) '  , dai_dqj_g ,  glob_cs_q.dai_dqj_cs);
        compare('(K8) '  , dIai_dqj_g ,  glob_cs_q.dIai_dqj_cs);
        compare('(K9) '  , dpsiddi_dqj_g ,  glob_cs_q.dpsiddi_dqj_cs);
        compare('(K10) '  , dBCi_dqj_g ,  glob_cs_q.dBCi_dqj_cs);
        compare('(K11) '  , dfi_dqj_g ,  glob_cs_q.dfi_dqj_cs);
        compare('(K12) '  , dfCi_dqj_g ,  glob_cs_q.dfCi_dqj_cs);
        compare('(K13) '  , dSiT_dqj_g , glob_cs_q.dSiT_dqj_cs);
        compare('(K14) '  , dSdi_dvj_g , glob_cs_v.dSdi_dvj_cs);
        compare('(K15) '  , dpsidi_dvj_g , glob_cs_v.dpsidi_dvj_cs);
        compare('(K16) '  , dBCi_dvj_g ,  glob_cs_v.dBCi_dvj_cs);
       

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
           error('%s is out of tolerance',txt);     
    else
        fprintf('%12s = %.3e  \x2713\n%s\n',txt,e);         
    end
    
end







