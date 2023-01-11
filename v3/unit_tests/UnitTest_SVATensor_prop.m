clear all; clc; 

%% Unit-tests for SVM Properties



Imat_vars = sym('I',[10 1],'real');
I = inertiaVecToMat(Imat_vars);                                                                        % Inertia matrix

u = sym('u',[6,1],'real');                                                                             % 6 x 1                                                                                                                                                                                                       % spatial-vector
v = sym('v',[6,1],'real');                                                                             % 6 x 1
w = sym('w',[6,1],'real');                                                                             % 6 x 2 
f = sym('f',[6,1],'real');                                                                             % 6 x 2

                                                                                                       % spatial matrices
U = sym('U',[6,randi([2,6])],'real');                                                                  % 6 x n1
V = sym('V',[6,randi([2,6])],'real');                                                                  % 6 x n2
W = sym('W',[6,randi([2,6])],'real');                                                                  % 6 x n3
F = sym('F',[6,randi([2,6])],'real');                                                                  % 6 x n4

                                                                                                       % scalar
syms lambda
n1 = randi([2,6]); n2 = randi([2,6]); n3 = randi([2,6]); n4 = randi([2,6]);
A = sym('U',[n1,n2],'real');                                                                           % n1 x n2 matrix
Y = sym('U',[n2,n3,n4],'real');                                                                        % n2 x n3 x n4 tensor

M1=simplify( cmfM(U)-(-rotT(crmM(U))));                                                                % M1
M2=simplify(-mT(V.',cmfM(U))-(rotT(Tm(crmM(U),V))));                                                   % M2
M3=simplify(Tm(-mT(V.',cmfM(U)),F)-(Tm(rotT(Tm(crmM(U),V)),F)));                                       % M3
M4=simplify(rotR(Tm(crmM(U),v))-(-crm(v)*U));                                                          % M4
M5=simplify( Tm(cmfM(U),F)-(rotR(Tm(cmf_barM(F),U))));                                                 % M5

M6=simplify(Tm(cmf_barM(F),U)-(rotR(Tm((cmfM(U)),F))));                                                % M6
M7=simplify(crmM(lambda*U)-(lambda*crmM(U)));                                                          % M7
M8=simplify(Tm(crmM(U),V)-(-rotR(Tm(crmM(V),U))));                                                     % M8
M9=simplify(crmM(crm(v)*U)-( mT(crm(v),crmM(U)) - Tm(crmM(U),crm(v))));                                % M9
M10=simplify(cmfM(crm(v)*U)-( mT(crf(v),cmfM(U)) - Tm(cmfM(U),crf(v))));                               % M10

M11=simplify(cmf_barM(rotR(Tm(cmfM(U),v)))-(Tm(cmfM(U),cmf_bar(v)) - mT(cmf_bar(v),crmM(U))));         % M11
M12=simplify(rotT(Tm(cmfM(U),F))-(-mT(F.',crmM(U))));                                                  % M12
M13_1=simplify( mT(V.',Tm(cmfM(U),F))-( Tm(rotT(rotR(Tm(crmM(V),U))),F)));                             % M13-1
M13_2=simplify(mT(V.',Tm(cmfM(U),F))-(  rotT(mT(F.',rotR(Tm(crmM(V),U)))) ));                          % M13-2
M14 =simplify( crf(v)*F-rotR( Tm(cmf_barM(F),v) ));                                                    % M14
M15 =simplify(cmf_bar(f)*U-rotR( Tm(cmfM(U),f) ));                                                     % M15

M16 =simplify(mT(V.',rotR(Tm(cmfM(U),F)))-(rotR(Tm(rotT(rotR(Tm(crmM(V),U))),F)) ));                   % M16
M17 =simplify(mT(V.',rotR(Tm(cmfM(U),F)))-(-rotT(mT(U.',rotR(Tm(cmfM(V),F)))) ));                      % M17
M18 =simplify(mT(I,rotR(Tm(cmfM(U),F)))-( rotR(mT(I,Tm(cmfM(U),F))) ));                                % M18
M19 = simplify(mT(I,rotR(Tm(crmM(U),V)))-( rotR(mT(I,Tm(crmM(U),V))) ));                               % M19
M20 =simplify(rotT(mT_v2(A,Y))-( Tm_v4(rotT(Y),A.') ));                                                % M20
M21 =simplify(mT_v2(F.',Tm_v4(crmM(U),V))-( -rotT(mT_v2(V.',rotR(Tm_v4(cmf_barM(F),U)))) ) );          % M21

M22 =simplify(Bmat(I,v).'*w-( -Bmat(I,w).'*v ) );                                                      % M22                                             
M23 = simplify(Bmat(I,v)*w - (Bmat(I,w)*v - I*crm(v)*w ) );                                            % M23
M24 = simplify( u.'*(Bmat(I,v)*w) - (-(v.'*Bmat(I,u)*w).') );                                          % M24

M25 =simplify( Tm_v4(rotT(Bten(I,V)),W)   -( -rotR(Tm_v4(rotT(Bten(I,W)),V)) ) );                      % M25
M26 =simplify( Tm_v4(Bten(I,V),W)   -(rotR(Tm_v4(Bten(I,W),V))-mT_v2(I,Tm_v4(crmM(V),W)) ) );          % M26
M27 =simplify( mT_v2(U.',rotR(Tm_v4(Bten(I,V),W)))   -(-rotT(mT_v2(V.',rotR(Tm_v4(Bten(I,U),W)))) ) ); % M27


compare('(M1) '  , double(~has(M1,sym(0))) );
compare('(M2) '  , double(~has(M2,sym(0))) );
compare('(M3) '  , double(~has(M3,sym(0))) );
compare('(M4) '  , double(~has(M4,sym(0))) );
compare('(M5) '  , double(~has(M5,sym(0))) );
compare('(M6) '  , double(~has(M6,sym(0))) );
compare('(M7) '  , double(~has(M7,sym(0))) );
compare('(M8) '  , double(~has(M8,sym(0))) );
compare('(M9) '  , double(~has(M9,sym(0))) );
compare('(M10) '  , double(~has(M10,sym(0))) );
compare('(M11) '  , double(~has(M11,sym(0))) );
compare('(M12) '  , double(~has(M12,sym(0))) );
compare('(M13 : 1) '  , double(~has(M13_1,sym(0))) );
compare('(M13 : 2) '  , double(~has(M13_2,sym(0))) );
compare('(M14) '  , double(~has(M14,sym(0))) );
compare('(M15) '  , double(~has(M15,sym(0))) );
compare('(M16) '  , double(~has(M16,sym(0))) );
compare('(M17) '  , double(~has(M17,sym(0))) );
compare('(M18) '  , double(~has(M18,sym(0))) );
compare('(M19) '  , double(~has(M19,sym(0))) );
compare('(M20) '  , double(~has(M20,sym(0))) );
compare('(M21) '  , double(~has(M21,sym(0))) );
compare('(M22) '  , double(~has(M22,sym(0))) );
compare('(M23) '  , double(~has(M23,sym(0))) );
compare('(M24) '  , double(~has(M24,sym(0))) );
compare('(M25) '  , double(~has(M25,sym(0))) );
compare('(M26) '  , double(~has(M26,sym(0))) );
compare('(M27) '  , double(~has(M27,sym(0))) );

%% Functions
function compare(txt, v1)
    if (size(v1,3))>1 
        e = tens_norm(v1);
    else
        e = norm(v1);
    end
    if e > 1e-10
        x = 'X';
        fprintf('%12s = %.3e  %s\n',txt,e,x);     
    else
        fprintf('%12s = %.3e  \x2713\n%s\n',txt,e);         
    end
    
end