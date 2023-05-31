% Unit tests for ID SO expressions

clear
N = 5;

% Create a random model with N links
model = autoTree(N, 1.5, pi/3);

checkDerivatives(model,'Floating Base No Rotors');


function checkDerivatives(model, desc)
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
m6 = rand(6,1); % random 6-motion vector
f6 = rand(6,1); %random 6-force vector

%% complex-step

newConfig = @(x) configurationAddition(model,q,x);
[out_cs.ddtau2_q_cs, out_cs.ddtau_dvq_cs] = ...
        complexStepHessian(@(x) ID_derivatives(model, newConfig(x) ,qd ,qdd),zeros(model.NV,1));
[~, out_cs.ddtau2_v_cs] = ...
        complexStepHessian(@(x) ID_derivatives(model, q ,x ,qdd) , qd);
out_cs.M_FO = complexStepJacobian_MFO(@(x) CRBA(model,  newConfig(x)), zeros(model.NV,1));

[glob,bod] = GlobalDynamics( model, q, qd, qdd,m6,f6);

fprintf('Running all possible cases for i,j,k \n')

for ii=1:model.NB
    for jj=1:ii
        for kk=1:jj
        
        ni = model.nv(ii); nj = model.nv(jj); nk = model.nv(kk);
        out_cs.i_ind = model.vinds{ii};   out_cs.j_ind = model.vinds{jj};   out_cs.k_ind = model.vinds{kk}; 

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
        
        A1 = Tm(cmfM(Si),BCi)-mT(BCi,crmM(Si));
        A2 = Tm(cmfM(Si),ICi)-mT(ICi,crmM(Si));
       
        if (ismember(kk,model.ancestors{jj})&&ismember(jj,model.ancestors{ii}))                                % k<=j<=i
           fprintf("\n ii = % d; jj= %d; kk= %d  \n \n",ii,jj,kk)              
 
            out_ana.SO_q_1 = -rotT( mT(psidj.',rotR(Tm(2*BIic_Si,psidk)))+...
                     mT(2*Sj.',rotR(Tm(cmf_barM(BCi.'*Si),psidk)))+...
                     mT(Sj.',rotR(Tm(cmf_barM(ICi*Si),psiddk))));                                              % SO q -1
            out_ana.SO_vq_1 = -rotT(mT(Sj.',rotR(Tm(2*BIic_Si,psidk))));                                       % SO vq -1   
            out_ana.MFO_1 = zeros(nj,ni,nk);                                                                   % M FO -1  
            out_ana.MFO_2 = zeros(ni,nj,nk);                                                                   % M FO -2 
            
            run_comp({'q-1','vq-1','MFO-1','MFO-2'},out_cs,out_ana);
                 
           if kk~=jj                                                                                           % k<j<=i
            out_ana.SO_q_2 = rotR(out_ana.SO_q_1);                                                             % SO q -2
            out_ana.SO_q_3 = mT(Sk.',rotR(Tm(2*BIic_psidi+2*A1,psidj)+Tm(A2,psiddj))+...
                      Tm(cmfM(Sj),2*BCi*psidi+ICi*psiddi+cmf_bar(fCi)*Si));                                    % SO q -3
            out_ana.SO_v_1 = -rotT(mT(Sj.',rotR(Tm(2*BIic_Si,Sk))));                                           % SO v -1
            out_ana.SO_v_2 = rotR(out_ana.SO_v_1);                                                             % SO v -2
            out_ana.SO_vq_3 = rotT(mT(Sk.',rotR(Tm(2*(-BIic_Si+cmf_barM(ICi*Si)),psidj)+2*Tm(cmf_barM(BCi.'*Si),Sj) ))+...
                            mT((psidk+Sdk).',rotR(Tm(cmf_barM(ICi*Si),Sj))));                                  % SO vq -3   
            out_ana.SO_vq_4 = mT(Sk.',rotR(Tm(2*BIic_Si,psidj)+Tm(cmf_barM(2*BCi*Si+ICi*(psidi+Sdi)),Sj ) ) ); % SO vq -4   
            out_ana.MFO_3 = mT(Sk.',rotR(Tm(cmf_barM(ICi*Si),Sj) ));                                           % M FO -3  
            out_ana.MFO_4 = rotT(out_ana.MFO_3);                                                               % M FO -4  

             run_comp({'q-2','q-3','v-1','v-2','vq-3','vq-4','MFO-3','MFO-4'},out_cs,out_ana);
               
             if jj~=ii                                                                                         % k<j<i
               out_ana.SO_q_4 = rotR(out_ana.SO_q_3);                                                          % SO q -4
               out_ana.SO_v_4 = mT(Sk.',rotR(Tm(2*BIic_Si,Sj)));                                               % SO v -4
               out_ana.SO_v_5 = rotR(out_ana.SO_v_4);                                                          % SO v -5
               out_ana.SO_vq_6 = mT(Sk.',Tm(2*(BIic_psidi+A1),Sj)+Tm(A2,psidj+Sdj));                           % SO vq -6  
               
               run_comp({'q-4','v-4','v-5','vq-6'},out_cs,out_ana);
               
             else                                                                                              % k<j=i
              out_ana.SO_v_6 = mT(Sk.',rotR( Tm(Tm(cmfM(Si),ICi)+cmf_barM(ICi*Si),Sj) ) );                     % SO v -6
              run_comp({'v-6'},out_cs,out_ana);
       
             end
           else                                                                                                % k=j<=i
              out_ana.SO_v_3 = -rotT(mT(Sj.',rotR(Tm(A2,Sk))));                                                % SO v -3
              run_comp({'v-3'},out_cs,out_ana);
           end
           
           if jj~=ii                                                                                           % k<=j<i
            out_ana.SO_q_5 = mT(Sj.',2*Tm(BIic_psidi+A1,psidk)+ Tm(A2,psiddk));                                % SO q -5
            out_ana.SO_q_6 = rotR(out_ana.SO_q_5);                                                             % SO q -6
            out_ana.SO_v_7 = mT(Sj.',Tm(2*BIic_Si,Sk));                                                        % SO v -7
            out_ana.SO_v_8 = rotR(out_ana.SO_v_7);                                                             % SO v -8
            out_ana.SO_vq_2 = mT(Sj.',rotR(Tm(2*BIic_Si,psidk)));                                              % SO vq -2
            out_ana.SO_vq_5 = mT(Sj.',Tm(2*(BIic_psidi+A1),Sk)+Tm(A2,psidk+Sdk));                              % SO vq -5
            out_ana.MFO_5 = mT(Sk.',Tm(A2,Sj) );                                                               % M FO -5  
            out_ana.MFO_6 = rotT(out_ana.MFO_5);                                                               % M FO -6  
          
            run_comp({'q-5','q-6','v-7','v-8','vq-2','vq-5','MFO-5','MFO-6'},out_cs,out_ana);
           end
        else
           fprintf("\n ii =% d; jj= %d; kk= %d  \n \n",ii,jj,kk)              
           fprintf("No Case Running\n")              

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
    else
        fprintf('%12s = %.3e  \x2713\n%s\n',txt,e);         
    end
    
end

function run_comp(case_str,out_cs,out_ana)
    
      for ii=1:numel(case_str)
          
        switch case_str{ii} 
         %q-cases   
         case 'q-1'
         compare('(SO (q) - 1)'  , out_ana.SO_q_1 , out_cs.ddtau2_q_cs(out_cs.i_ind,out_cs.j_ind,out_cs.k_ind));
         case 'q-2'
         compare('(SO (q) - 2)'  , out_ana.SO_q_2 , out_cs.ddtau2_q_cs(out_cs.i_ind,out_cs.k_ind,out_cs.j_ind));
         case 'q-3'                     
         compare('(SO (q) - 3)'  , out_ana.SO_q_3 , out_cs.ddtau2_q_cs(out_cs.k_ind,out_cs.i_ind,out_cs.j_ind));
         case 'q-4'                            
         compare('(SO (q) - 4)'  , out_ana.SO_q_4 , out_cs.ddtau2_q_cs(out_cs.k_ind,out_cs.j_ind,out_cs.i_ind));
         case 'q-5'        
         compare('(SO (q) - 5)'  , out_ana.SO_q_5 , out_cs.ddtau2_q_cs(out_cs.j_ind,out_cs.k_ind,out_cs.i_ind));
         case 'q-6'             
         compare('(SO (q) - 6)'  , out_ana.SO_q_6 , out_cs.ddtau2_q_cs(out_cs.j_ind,out_cs.i_ind,out_cs.k_ind));
          
        %v-cases   
        case 'v-1'                            
        compare('(SO (v) - 1)'  , out_ana.SO_v_1 , out_cs.ddtau2_v_cs(out_cs.i_ind,out_cs.j_ind,out_cs.k_ind));
        case 'v-2'                            
        compare('(SO (v) - 2)'  , out_ana.SO_v_2 , out_cs.ddtau2_v_cs(out_cs.i_ind,out_cs.k_ind,out_cs.j_ind));
        case 'v-3'                            
        compare('(SO (v) - 3)'  , out_ana.SO_v_3 , out_cs.ddtau2_v_cs(out_cs.i_ind,out_cs.j_ind,out_cs.k_ind));
        case 'v-4'                            
        compare('(SO (v) - 4)'  , out_ana.SO_v_4 , out_cs.ddtau2_v_cs(out_cs.k_ind,out_cs.i_ind,out_cs.j_ind));
        case 'v-5'                            
        compare('(SO (v) - 5)'  , out_ana.SO_v_5 , out_cs.ddtau2_v_cs(out_cs.k_ind,out_cs.j_ind,out_cs.i_ind));
        case 'v-6'                            
        compare('(SO (v) - 6)'  , out_ana.SO_v_6 , out_cs.ddtau2_v_cs(out_cs.k_ind,out_cs.i_ind,out_cs.j_ind));
        case 'v-7'                            
        compare('(SO (v) - 7)'  , out_ana.SO_v_7 , out_cs.ddtau2_v_cs(out_cs.j_ind,out_cs.k_ind,out_cs.i_ind));
        case 'v-8'                            
        compare('(SO (v) - 8)'  , out_ana.SO_v_8 , out_cs.ddtau2_v_cs(out_cs.j_ind,out_cs.i_ind,out_cs.k_ind));
                 
        %qv-cases   
        case 'vq-1'                            
        compare('(SO(vq) - 1)'  , out_ana.SO_vq_1 , out_cs.ddtau_dvq_cs(out_cs.i_ind,out_cs.j_ind,out_cs.k_ind));
        case 'vq-2'                            
        compare('(SO(vq) - 2)'  , out_ana.SO_vq_2 , out_cs.ddtau_dvq_cs(out_cs.j_ind,out_cs.i_ind,out_cs.k_ind));
        case 'vq-3'                            
        compare('(SO(vq) - 3)'  , out_ana.SO_vq_3 , out_cs.ddtau_dvq_cs(out_cs.i_ind,out_cs.k_ind,out_cs.j_ind));
        case 'vq-4'                            
        compare('(SO(vq) - 4)'  , out_ana.SO_vq_4 , out_cs.ddtau_dvq_cs(out_cs.k_ind,out_cs.i_ind,out_cs.j_ind));
        case 'vq-5'                            
        compare('(SO(vq) - 5)'  , out_ana.SO_vq_5 , out_cs.ddtau_dvq_cs(out_cs.j_ind,out_cs.k_ind,out_cs.i_ind));
        case 'vq-6'                            
        compare('(SO(vq) - 6)'  , out_ana.SO_vq_6 , out_cs.ddtau_dvq_cs(out_cs.k_ind,out_cs.j_ind,out_cs.i_ind));
         
        % aq Cases
        case 'MFO-1'                            
        compare('(FO(M) - 1)'  , out_ana.MFO_1 , out_cs.M_FO(out_cs.j_ind,out_cs.i_ind,out_cs.k_ind));
        case 'MFO-2'                            
        compare('(FO(M) - 2)'  , out_ana.MFO_2 , out_cs.M_FO(out_cs.i_ind,out_cs.j_ind,out_cs.k_ind));
        case 'MFO-3'                            
        compare('(FO(M) - 3)'  , out_ana.MFO_3 , out_cs.M_FO(out_cs.k_ind,out_cs.i_ind,out_cs.j_ind));
        case 'MFO-4'                            
        compare('(FO(M) - 4)'  , out_ana.MFO_4 , out_cs.M_FO(out_cs.i_ind,out_cs.k_ind,out_cs.j_ind));
        case 'MFO-5'                            
        compare('(FO(M) - 5)'  , out_ana.MFO_5 , out_cs.M_FO(out_cs.k_ind,out_cs.j_ind,out_cs.i_ind));
        case 'MFO-6'                            
        compare('FO(M) - 6)'  , out_ana.MFO_6 , out_cs.M_FO(out_cs.j_ind,out_cs.k_ind,out_cs.i_ind));
               
        
        otherwise
        fprintf('Error- No Case Selected')
        end
      end



end






