clc; clear all;

% run([pwd,'\..\startup.m'])
% Testing the FO partials of spatial force

N = 17;

% Create a random model with N links
model = autoTree(N, 3);
model.jtype{1}='Fb';
model = postProcessModel(model);

q  = ones(model.NQ,1);

if strcmp(model.jtype{1},'SE3')
    q(1:16) = reshape( randomSE3(),[16 1]);
elseif strcmp(model.jtype{1},'Fb')  % normalizing the quoternion
    q(1:4) = q(1:4)/norm(q(1:4));
end

qd = ones(model.NV,1);
qdd = ones(model.NV,1);

fprintf("N = %d\n", N)% 

fprintf("model.NQ = %d\n", model.NQ)% 
fprintf("model.NV = %d\n", model.NV)% 

%% complex-step- FO derivatives of spatial force/spatial cumulative force
step=1e-20; % step size for complex-step

newConfig = @(x) configurationAddition(model,q,x);

[glob,bod] = GlobalDynamics_red( model, q, qd, qdd);

stacked_f = zeros(6,model.NB);
stacked_df_dq = zeros(6,model.NB,model.NV);

stacked_fc = zeros(6,model.NB);
stacked_dfc_dq = zeros(6,model.NB,model.NV);


for ii=1:N
    for jj=1:N
        
        ni = model.nv(ii); nj = model.nv(jj); 
        ii_vec = model.vinds{ii}; jj_vec = model.vinds{jj};
           
          % Getting the partials of the global dynamics w.r.t q
         [glob_cs_q,bod_cs_q] = perturbIndex_q_sp(model,@(x) GlobalDynamics_red(model, newConfig(x) ...
            ,qd ,qdd), zeros(model.NV,1),jj,ii ); % complex with perturb q
       
          % Getting the partials of the global dynamics w.r.t v
          glob_cs_v = perturbIndex_v_sp(model,@(x) GlobalDynamics_red(model, q ...
            ,x ,qdd), qd,jj,ii );                 % complex with perturb v
        
          % Getting the partials of the global dynamics w.r.t a
          glob_cs_a = perturbIndex_a_sp(model,@(x) GlobalDynamics_red(model, q ...
            ,qd ,x), qd,jj,ii );                 % complex with perturb v
        
         % glob quantities
        Si_g = glob.S{ii};       psidi_g = glob.psid{ii};  vi_g = glob.v{ii};      ai_g = glob.a{ii};
        Sdi_g = glob.Sd{ii};     Ii_g = glob.I{ii};        Iai_g = glob.Ia{ii};    psiddi_g = glob.psidd{ii};
        vJxSi_g = glob.vJxS{ii}; ICi_g = glob.IC{ii};      BCi_g = glob.BC{ii};    fCi_g = glob.fC{ii}; 
        fi_g = glob.f{ii};        Bi_g = glob.B{ii}; 
        iXj =  glob.iXO{ii}*glob.OXi{jj};     
        iX0 = glob.iXO{ii};
        OXi = glob.OXi{ii};
        
        Sj_g = glob.S{jj};       psidj_g = glob.psid{jj};  psiddj_g = glob.psidd{jj};  
        Sdj_g = glob.Sd{jj};     ICj_g = glob.IC{jj};      BCj_g = glob.BC{jj};
        fCj_g = glob.fC{jj};
        BIic_psidj_g = Bten(ICi_g,psidj_g);                BIjc_psidj_g = Bten(ICj_g,psidj_g);
        BIic_Sj_g = Bten(ICi_g,Sj_g);       
        
        dfi_dqj_g = zeros(6,nj);      
        dfCi_dqj_g = zeros(6,nj);      
        stacked_f(:,ii) = fi_g;
        stacked_fc(:,ii) = fCi_g;

        if  ismember(jj,model.ancestors{ii})                                             % jj<=ii
         fprintf("\n ii =% d; jj= %d;  \n \n",ii,jj)
        
         dfi_dqj_g = Ii_g*psiddj_g+ cmf_barM(fi_g)*Sj_g + 2*Bi_g*psidj_g;                               dfi_dqj_b =  iX0*(dfi_dqj_g+crmM(fi_g)*Sj_g);      % K7  % K11
         stacked_df_dq(1:6,ii,jj_vec) = dfi_dqj_g;
         
         dfCi_dqj_g = ICi_g*psiddj_g + cmf_barM(fCi_g)*Sj_g + 2*BCi_g*psidj_g;                         % K12
         stacked_dfc_dq(1:6,ii,jj_vec) = dfCi_dqj_g;

        elseif ismember(ii,model.ancestors{jj}) 
          fprintf("\n ii = %d; jj= %d  j>i Case \n \n",ii,jj)  
          
          dfCi_dqj_g = ICj_g*psiddj_g + cmf_barM(fCj_g)*Sj_g + 2*BCj_g*psidj_g;                        % K12
          stacked_dfc_dq(1:6,ii,jj_vec) = dfCi_dqj_g;

       else
         fprintf("\n ii = %d; jj= %d  different branches  \n \n",ii,jj)  

        end
        
        compare('(K11) '  , dfi_dqj_g ,  glob_cs_q.dfi_dqj_cs);
        compare('(K12) '  , dfCi_dqj_g ,  glob_cs_q.dfCi_dqj_cs);
        
        % Analytical algorithm to get the spatial force FO partials

        [df_dq_ana, dfc_dq_ana,df_dv_ana,dfc_dv_ana,df_da_ana,dfc_da_ana] = ...
                        spatial_force_derivatives( model, q, qd, qdd, ii,jj );
        
        compare('(df_dq ana) '  , df_dq_ana ,  glob_cs_q.dfi_dqj_cs);
        compare('(dfc_dq ana) '  , dfc_dq_ana ,  glob_cs_q.dfCi_dqj_cs);
        
        compare('(df_dv ana) '  , df_dv_ana ,  glob_cs_v.dfi_dqv_cs);
        compare('(dfc_dv ana) '  , dfc_dv_ana ,  glob_cs_v.dfci_dqv_cs);
        
        compare('(df_da ana) '  , df_da_ana ,  glob_cs_a.dfi_da_cs);
        compare('(dfc_da ana) '  , dfc_da_ana ,  glob_cs_a.dfci_da_cs);

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
    else
        fprintf('%12s = %.3e  \x2713\n%s\n',txt,e);         
    end
    
  %  fprintf('%12s = %.3e  %s\n',txt,e,x); 
end
    
    
    
    
    