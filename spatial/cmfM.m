function tens = cmfM(F)

% Description- 
% cmfM operator on a spatial force matrix/vector
%(F) x~*
% Inputs 
% F (6XN spatial force matrix)
% Outputs
% 6x6xN tensor

if (size(F,2)>1) % if F is a matrix
    
for ii=1:size(F,2)
      t1 = crf(F(:,ii));
      tens(:,:,ii) = t1;
end
 else % if F is a vector
   tens = crf(F);
end