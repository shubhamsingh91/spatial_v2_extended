function tens = crmM(U)

% Description- 
% crmM operator on a spatial motion matrix/vector
%(U) x~
% Inputs 
% U (6XN spatial motion matrix)
% Outputs
% 6x6xN tensor

if size(U,2)>1 % if U is a matrix
    
    for ii=1:size(U,2)
         t1 = crm(U(:,ii));
          tens(:,:,ii) = t1;
    end
else  % if U is a vector
     tens = crm(U);
end
    



