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
%$Id: friction.m 79 2017-04-20 11:07:54Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/friction.m $
%
%friction is a function that computed the dimensionless friction coefficient
%
%Cf=friction(h,Mak,Cf_old,input,fid_log,kt)
%
%INPUT:
%   -
%
%OUTPUT:
%   -
%
%HISTORY:
%160223
%   -V. Created for the first time.

function Cf=friction(h,Mak,Cf_old,input,fid_log,kt)
%comment out fot improved performance if the version is clear from github
% version='1';
% if kt==1; fprintf(fid_log,'friction version: %s\n',version); end 


%%

switch input.mdv.frictiontype
    case 1 %constant friction
        Cf=Cf_old;
    case 2 %friction related to grain size
        error('not yet implemented')
    case 3 %friction related to flow depth
        error('not yet implemented')        
end
