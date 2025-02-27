% Unit tests for multi-DoF Identities

clear
N = 5;

% Create a random model with N links
model = autoTree(N, 1.5, pi/3);
checkSVATensor_iden(model,'Floating Base No Rotors');

function checkSVATensor_iden(model, desc)
fprintf('====================================\n');
fprintf('%s\n',desc);
fprintf('====================================\n');

model.jtype{1}='SE3';
model = postProcessModel(model);

q  = rand(model.NQ,1);

if strcmp(model.jtype{1},'SE3')
    q(1:16) = reshape( randomSE3(),[16 1]);
elseif strcmp(model.jtype{1},'Fb')  % normalizing the quoternion
    q(1:4) = q(1:4)/norm(q(1:4));
end

qd = rand(model.NV,1);
qdd=rand(model.NV,1);

%%
step=1e-20; % step size for complex-step

newConfig = @(x) configurationAddition(model,q,x);

m6 = rand(6,1);                                                              % random 6-motion vector
f6 = rand(6,1);                                                              % random 6-force vector


[glob,bod] = GlobalDynamics( model, q, qd, qdd,m6,f6);
                                                                             % For all possible ii and jj cases

for ii=1:model.NB                           
     for jj=1:ii

        ni = model.nv(ii); nj = model.nv(jj);
        p = randi([1,nj]);                                                   % 1 random single mode of joint-j

        [glob_cs_q] = perturbIndex_q(model,@(x) GlobalDynamics(model, newConfig(x) ...
            ,qd ,qdd,m6,f6), zeros(model.NV,1),jj,ii,p );                    % complex with perturb q

        [glob_cs_v] = perturbIndex_v(model,@(x) GlobalDynamics(model, q ...
            ,x ,qdd,m6,f6), qd,jj,ii );                                      % complex with perturb v

                                                                             % glob coordinates
        Si_g = glob.S{ii};       psidi_g = glob.psid{ii};  vi_g = glob.v{ii};      ai_g = glob.a{ii};
        Sdi_g = glob.Sd{ii};     Ii_g = glob.I{ii};        Iai_g = glob.Ia{ii};    psiddi_g = glob.psidd{ii};
        vJxSi_g = glob.vJxS{ii}; ICi_g = glob.IC{ii};      BCi_g = glob.BC{ii};    fCi_g = glob.fC{ii}; 
        fi_g = glob.f{ii};        Bi_g = glob.B{ii}; 
        oXi_g = glob.OXi{ii};

        Sj_g = glob.S{jj};       psidj_g = glob.psid{jj};  psiddj_g = glob.psidd{jj};  
        Sdj_g = glob.Sd{jj};     ICj_g = glob.IC{jj};      BCj_g = glob.BC{jj};
        fCj_g = glob.fC{jj};
        BIic_psidj_g = Bten(ICi_g,psidj_g);                BIjc_psidj_g = Bten(ICj_g,psidj_g);
        BIic_Sj_g = Bten(ICi_g,Sj_g);                      BIjc_Sj_g = Bten(ICj_g,Sj_g);
        vj_g = glob.v{jj};
        vJj_g = glob.vJ{jj};
        vlamj_g = glob.vlam{jj};
        Sjp_g = glob.S{jj}(:,p);
        
        if model.parent(jj)==0
         xilamj_g = zeros(6,1);
         gammalamj_g = zeros(6,1);
        else
         xilamj_g = glob.xi{model.parent(jj)};
         gammalamj_g = glob.gamma{model.parent(jj)};
        end
        
        Ivi_g = glob.Iv{ii};
        Imi_g = glob.Im{ii};
        xii_g = glob.xi{ii};
        gammai_g = glob.gamma{ii};
        vixsf_g = glob.vxsf{ii};
        SiTf_g = glob.STf{ii};
         
        dIvi_dqj_g = zeros(6,nj);     
        dImi_dqj_g = zeros(6,nj);     
        dvisf_dqj_g = zeros(6,nj);     
        dSiTf_dqj_g = zeros(ni,nj);     
        dxii_dqj_g = zeros(6,nj);     
        dgammai_dqj_g = zeros(6,nj);   
        doXi_dqj_g = zeros(6,6,nj);
        
        dvi_dvj_g = zeros(6,nj);
        dxii_dvj_g = zeros(6,nj);
        dSi_dqjp_g = zeros(6,ni);
        
        if  ismember(jj,model.ancestors{ii})                                 % jj<=ii
            dSi_dqjp_g = crm(Sjp_g)*Si_g;                                    % J1
            dvisf_dqj_g = cmf_bar(f6)*crm(vj_g-vi_g-vJj_g)*Sj_g;             % J2
            dSiTf_dqj_g = -Si_g.'*cmf_bar(f6)*Sj_g;                          % J3
            dImi_dqj_g = cmf_bar(Imi_g)*Sj_g + Ii_g*crm(m6)*Sj_g;            % J4
            dIvi_dqj_g = cmf_bar(Ivi_g)*Sj_g + Ii_g*psidj_g;                 % J5
            dxii_dqj_g = crm(vlamj_g-vi_g)*psidj_g+crm(xilamj_g-xii_g)*Sj_g; % J6
            dgammai_dqj_g = crm(gammalamj_g-gammai_g)*Sj_g;                  % J7
            dvi_dvj_g = Sj_g;                                                % J8
            dxii_dvj_g = crm(vlamj_g-vi_g)*Sj_g + Sdj_g;                     % J9
          
            if jj==ii
            fprintf("\n ii = %d; jj= %d  j<=i Case \n \n",ii,jj)
            else
              fprintf("\n ii = %d; jj= %d  j<i Case \n \n",ii,jj) 
            end

            elseif ismember(ii,model.ancestors{jj}) 
                 fprintf("\n ii = % d; jj= %d  j>i Case \n \n",ii,jj)  
               
            else
                  fprintf("\n ii = % d; jj= %d  different branches  \n \n",ii,jj)  

        end
        fprintf('Global Coordinates\n');
        
        compare('(J1) '  , dSi_dqjp_g ,  glob_cs_q.dSi_dqjp_cs);    
        compare('(J2) '  , dvisf_dqj_g ,  glob_cs_q.dvisf_dqj_cs);    
        compare('(J3) '  , dSiTf_dqj_g ,  glob_cs_q.dSiTf_dqj_cs);
        compare('(J4) '  , dImi_dqj_g ,  glob_cs_q.dImi_dqj_cs);    
        compare('(J5) '  , dIvi_dqj_g ,  glob_cs_q.dIvi_dqj_cs);
        compare('(J6) '  , dxii_dqj_g ,  glob_cs_q.dxii_dqj_cs);    
        compare('(J7) '  , dgammai_dqj_g ,  glob_cs_q.dgammai_dqj_cs);             
        compare('(J8) '  , dvi_dvj_g ,  glob_cs_v.dvi_dvj_cs);             
        compare('(J9) '  , dxii_dvj_g ,  glob_cs_v.dxii_dvj_cs);             

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







