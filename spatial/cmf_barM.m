function tens = cmf_barM(F)

% Description- 
% cmf_bar operator on a spatial force matrix/vector
%(F) x~-*
% Inputs 
% F (6XN spatial force matrix)
% Outputs
% 6x6xN tensor

if (size(F,2)>1)
for ii=1:size(F,2)
  tens(:,:,ii) =  cmf_bar(F(:,ii));
end
else
   tens = cmf_bar(F); 
end

