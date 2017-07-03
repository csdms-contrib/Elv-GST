%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 7 $
%$Date: 2017-02-06 14:42:31 +0100 (Mon, 06 Feb 2017) $
%$Author: V $
%$Id: exist_str.m 7 2017-02-06 13:42:31Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/auxiliary/exist_str.m $
%
function isFieldResult = exist_str (inStruct, fieldName)
% inStruct is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
isFieldResult = 0;
f = fieldnames(inStruct(1));
for i=1:length(f)
if(strcmp(f{i},strtrim(fieldName)))
isFieldResult = 1;
return;
elseif isstruct(inStruct(1).(f{i}))
isFieldResult = exist_str(inStruct(1).(f{i}), fieldName);
if isFieldResult
return;
end
end
end