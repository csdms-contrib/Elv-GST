%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 79 $
%$Date: 2017-04-20 13:07:54 +0200 (Thu, 20 Apr 2017) $
%$Author: V $
%$Id: ini_La.m 79 2017-04-20 11:07:54Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/ini_La.m $
%
%ini_La this function creates the initial active layer thickness vector
%
%La=ini_La(input,fid_log)
%
%INPUT:
%   -input = variable containing the input [struct] e.g. input
%
%OUTPUT:
%   -
%
%HISTORY:
%160223
%   -V. Created for the first time.

function La=ini_La(input,fid_log)
%comment out fot improved performance if the version is clear from github
% version='1';
% fprintf(fid_log,'ini_La version: %s\n',version);

switch input.mor.Latype
    case {1,4} %constant active layer thickness
        La=repmat(input.mor.La,1,input.mdv.nx);               
end
