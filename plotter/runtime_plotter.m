clc; clear all; 
close all;
fclose('all')

% For plotting TRO run-time results

%% Bar chart for pinocchio Timings

multi_dof_str_all = {'double_pendulum','ur3_robot','hyq_f','atlas_f','talos_full_v2_f'};
multi_dof_str_fixed = {'double_pendulum','ur3_robot','N','N','N'};

n_vec = [1:1: size(multi_dof_str_all,2)];

 
%% WINDOWS PATH in WINDOWS
% digits 16
% format longg
% GCC
 folderstr_gcc_avx_AD_v1 =  'data\AD_v1\avx\gcc\mdof\';
 folderstr_gcc_avx_AD_v2 =  'data\AD_v2\avx\gcc\mdof\';

 folderstr_gcc_avx_AD_v3 =  'data\AD_v3\avx\gcc\mdof\';
 folderstr_gcc_avx_AD_v4 =  'data\AD_v4\avx\gcc\mdof\';

 % CLANG
 folderstr_c_avx_AD_v1 =  'data\AD_v1\avx\clang\mdof\';
 folderstr_c_avx_AD_v2 =   'data\AD_v2\avx\clang\mdof\';
 
 folderstr_c_avx_AD_v3 =  'data\AD_v3\avx\clang\mdof\';
 folderstr_c_avx_AD_v4 =  'data\AD_v4\avx\clang\mdof\';
 
 % archive data
  folderstr_gcc_avx_AD_v1_arch =  'data\AD_v1\avx\gcc\mdof\';
 folderstr_gcc_avx_AD_v2_arch =  'data\AD_v2\avx\gcc\mdof\';

 folderstr_gcc_avx_AD_v3_arch =  'data\AD_v3\avx\gcc\mdof\';
 folderstr_gcc_avx_AD_v4_arch =  'data\AD_v4\avx\gcc\mdof\';

 % CLANG
 folderstr_c_avx_AD_v1_arch =  'data\AD_v1\avx\clang\mdof\';
 folderstr_c_avx_AD_v2_arch =   'data\AD_v2\avx\clang\mdof\';
 
 folderstr_c_avx_AD_v3_arch =  'data\AD_v3\avx\clang\mdof\';
 folderstr_c_avx_AD_v4_arch =  'data\AD_v4\avx\clang\mdof\';
 
%%
%%

 

%%
for ii=1:numel(multi_dof_str_all)
   %% GCC
       
    % v1
    str1 = multi_dof_str_all(ii); 
    str2 = strcat(folderstr_gcc_avx_AD_v1,str1,'.txt');
    fileID = fopen(str2{1});
    num_in = fscanf(fileID,'%f');       
          
    
    time_RNEA_deriv_FO_faster_avx_g(ii) = num_in(1);          % RNEA FO deriv based on IDSVA FO
    time_RNEA_deriv_FO_AD_cppad_nocg_avx_g(ii) = num_in(2);   % ID FO partials using cppad w/out cg
    time_RNEA_deriv_FO_AD_cas_nocg_avx_g(ii) = num_in(3);     % ID FO partials using casadi w/out cg
    time_RNEA_deriv_SO_avx_g(ii) = num_in(4);                 % IDSVA SO analytical
    time_RNEA_deriv_SO_AD_cppad_nocg_avx_g(ii) = num_in(5);   % ID SO partials using cppad w/out cg
    time_RNEA_deriv_SO_AD_cas_nocg_avx_g(ii) = num_in(6);     % ID SO partials using casadi w/out cg

     
    %v2
    str2_v2 = strcat(folderstr_gcc_avx_AD_v2,str1,'.txt');
    fileID_v2 = fopen(str2_v2{1});
    num_in_v2 = fscanf(fileID_v2,'%f');  
    
    time_RNEA_deriv_FO_AD_cas_cg_avx_g(ii) = num_in_v2(2);   % ID FO partials using casadi w cg
    time_RNEA_deriv_SO_AD_cas_cg_avx_g(ii) = num_in_v2(3);   % ID SO partials using casadi w cg
    time_RNEA_deriv_SO_ana_cg_avx_g(ii) =    num_in_v2(4);   % ID SO partials using IDSVA w cg

     %v2- arch
    str2_v2 = strcat(folderstr_gcc_avx_AD_v2_arch,str1,'.txt');
    fileID_v2 = fopen(str2_v2{1});
    num_in_v2_arch = fscanf(fileID_v2,'%f');  
    
    time_RNEA_deriv_FO_AD_cas_cg_avx_g_arch(ii) = num_in_v2_arch(2);   % ID FO partials using casadi w cg
    time_RNEA_deriv_SO_AD_cas_cg_avx_g_arch(ii) = num_in_v2_arch(3);   % ID SO partials using casadi w cg
    time_RNEA_deriv_SO_ana_cg_avx_g_arch(ii) =    num_in_v2_arch(4);   % ID SO partials using IDSVA w cg
    
    % v3
    str2_v3 = strcat(folderstr_gcc_avx_AD_v3,str1,'.txt');
    fileID_v3 = fopen(str2_v3{1});
    num_in_v3 = fscanf(fileID_v3,'%f');  
   
    time_ABA_deriv_FO_avx_g(ii) = num_in_v3(1);               % ABA FO deriv based on FDSVA
  

    % v4
    str2_v4 = strcat(folderstr_gcc_avx_AD_v4,str1,'.txt');
    fileID_v4 = fopen(str2_v4{1});
    num_in_v4 = fscanf(fileID_v4,'%f'); 
    
    time_ABA_deriv_FO_AD_cas_cg_avx_g(ii)   = num_in_v4(2); % FD FO partials using casadi w cg
    time_ABA_deriv_SO_avx_g(ii)             = num_in_v4(3); % FD SO deriv based on SVA
    time_ABA_deriv_SO_AD_cas_cg_avx_g(ii)   = num_in_v4(4); % FD SO partials using casadi w cg
    time_ABA_deriv_SO_ana_cg_avx_g(ii)      = num_in_v4(5); % FD SO partials using Analytical w cg

 
    %% CLANG
       
    % v1
    str1c = multi_dof_str_all(ii); 
    str2c = strcat(folderstr_c_avx_AD_v1,str1c,'.txt');
    fileIDc = fopen(str2c{1});
    num_inc = fscanf(fileIDc,'%f');       
          
    
    time_RNEA_deriv_FO_faster_avx_c(ii) = num_inc(1);          % RNEA FO deriv based on IDSVA FO
    time_RNEA_deriv_FO_AD_cppad_nocg_avx_c(ii) = num_inc(2);   % ID FO partials using cppad w/out cg
    time_RNEA_deriv_FO_AD_cas_nocg_avx_c(ii) = num_inc(3);     % ID FO partials using casadi w/out cg
    time_RNEA_deriv_SO_avx_c(ii) = num_inc(4);                 % IDSVA SO analytical
    time_RNEA_deriv_SO_AD_cppad_nocg_avx_c(ii) = num_inc(5);   % ID SO partials using cppad w/out cg
    time_RNEA_deriv_SO_AD_cas_nocg_avx_c(ii) = num_inc(6);     % ID SO partials using casadi w/out cg

    
    %v2
    str2_v2c = strcat(folderstr_c_avx_AD_v2,str1,'.txt');
    fileID_v2c = fopen(str2_v2c{1});
    num_in_v2c = fscanf(fileID_v2c,'%f');  
    
    time_RNEA_deriv_FO_AD_cas_cg_avx_c(ii) = num_in_v2c(2);   % ID FO partials using casadi w cg
    time_RNEA_deriv_SO_AD_cas_cg_avx_c(ii) = num_in_v2c(3);   % ID SO partials using casadi w cg
    time_RNEA_deriv_SO_ana_cg_avx_c(ii) =    num_in_v2c(4);   % ID SO partials using IDSVA w cg


%     % v3
    str2_v3c = strcat(folderstr_c_avx_AD_v3,str1,'.txt');
    fileID_v3c = fopen(str2_v3c{1});
    num_in_v3c = fscanf(fileID_v3c,'%f');  
   
    time_ABA_deriv_FO_avx_c(ii) = num_in_v3c(1);               % ABA FO deriv based on FDSVA
  

    
    %v4
    str2_v4c = strcat(folderstr_c_avx_AD_v4,str1,'.txt');
    fileID_v4c = fopen(str2_v4c{1});
    num_in_v4c = fscanf(fileID_v4c,'%f'); 
    
    time_ABA_deriv_FO_AD_cas_cg_avx_c(ii)   = num_in_v4c(2); % FD FO partials using casadi w cg
    time_ABA_deriv_SO_avx_c(ii)             = num_in_v4c(3); % FD SO deriv based on SVA
    time_ABA_deriv_SO_AD_cas_cg_avx_c(ii)   = num_in_v4c(4); % FD SO partials using casadi w cg
    time_ABA_deriv_SO_ana_cg_avx_c(ii)   = num_in_v4c(5); % FD SO partials using Analytical w cg


end

% 
n_full = numel(n_vec);
time_RNEA_deriv_FO_AD_cppad_cg_avx_g = zeros(1,n_full);
time_RNEA_deriv_SO_AD_cppad_cg_avx_g = zeros(1,n_full);
time_ABA_deriv_FO_AD_cppad_cg_avx_g = zeros(1,n_full);
time_ABA_deriv_SO_AD_cppad_cg_avx_g = zeros(1,n_full);

time_RNEA_deriv_FO_AD_cppad_cg_avx_c = zeros(1,n_full);
time_RNEA_deriv_SO_AD_cppad_cg_avx_c = zeros(1,n_full);
time_ABA_deriv_FO_AD_cppad_cg_avx_c = zeros(1,n_full);
time_ABA_deriv_SO_AD_cppad_cg_avx_c = zeros(1,n_full);
            
                        
%% ID SO results- using gcc- without no-codegen versions- TRO

close all
figure_pinocchio_bar_IDSVA_SO(time_RNEA_deriv_SO_avx_g,...
                    time_RNEA_deriv_SO_AD_cas_cg_avx_g,...
                    time_RNEA_deriv_SO_ana_cg_avx_g,...
                    time_RNEA_deriv_SO_AD_cas_nocg_avx_g,...
                    time_RNEA_deriv_SO_avx_c,...
                    time_RNEA_deriv_SO_AD_cas_cg_avx_c,...
                    time_RNEA_deriv_SO_ana_cg_avx_c,...
                    time_RNEA_deriv_SO_AD_cas_nocg_avx_c)    
                
   % IDSVA+CodeGen to IDSVA             
     speedups_IDSVA_SO_to_IDSVA_SO_cg_g = time_RNEA_deriv_SO_ana_cg_avx_g./time_RNEA_deriv_SO_avx_g;  % speedup of ID SO to ID SO cg
     speedups_IDSVA_SO_to_IDSVA_SO_cg_c = time_RNEA_deriv_SO_ana_cg_avx_c./time_RNEA_deriv_SO_avx_c;  % speedup of ID SO to ID SO cg

   % IDSVA to AD+cg            
      speedups_IDSVA_SO_to_AD_SO_cg_g = time_RNEA_deriv_SO_AD_cas_cg_avx_g./time_RNEA_deriv_SO_avx_g; % speedup of ID SO to ID SO cg
      speedups_IDSVA_SO_to_AD_SO_cg_c = time_RNEA_deriv_SO_AD_cas_cg_avx_c./time_RNEA_deriv_SO_avx_c;
      
   % AD+cg to AD           
     speedups_AD_cg_to_non_cg_g = time_RNEA_deriv_SO_AD_cas_nocg_avx_g./time_RNEA_deriv_SO_AD_cas_cg_avx_g;
     speedups_AD_cg_to_non_cg_c = time_RNEA_deriv_SO_AD_cas_nocg_avx_c./time_RNEA_deriv_SO_AD_cas_cg_avx_c;


 %% FD SO results- using gcc/clang- without no-codegen versions- for TRO
 
% FDSVA to AD+cg            
speedups_FDSVASO_to_FD_SO_AD_no_cg_g =time_ABA_deriv_SO_AD_cas_cg_avx_g./time_ABA_deriv_SO_avx_g;
speedups_FDSVASO_to_FD_SO_AD_no_cg_c =time_ABA_deriv_SO_AD_cas_cg_avx_c./time_ABA_deriv_SO_avx_c;

% FDSVA+cg to FDSVA            
% speedups_FDSVASO_to_FD_SO_AD_cg_g =time_ABA_deriv_SO_ana_cg_avx_g./time_ABA_deriv_SO_avx_g;
% speedups_FDSVASO_to_FD_SO_AD_cg_c =time_ABA_deriv_SO_ana_cg_avx_c./time_ABA_deriv_SO_avx_c;


% close all
figure_pinocchio_bar_FDSVA_SO(time_ABA_deriv_SO_avx_g,...
                    time_ABA_deriv_SO_AD_cas_cg_avx_g,...
                    time_ABA_deriv_SO_ana_cg_avx_g,...
                    time_ABA_deriv_SO_avx_c,...
                    time_ABA_deriv_SO_AD_cas_cg_avx_c,...
                    time_ABA_deriv_SO_ana_cg_avx_c)
                               
    
                
                