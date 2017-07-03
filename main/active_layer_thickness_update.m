%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 106 $
%$Date: 2017-06-27 14:45:13 +0200 (Tue, 27 Jun 2017) $
%$Author: V $
%$Id: active_layer_thickness_update.m 106 2017-06-27 12:45:13Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/active_layer_thickness_update.m $
%
%active_layer_thickness_update is a function that updates the active layer thickness
%
%\texttt{La=active_layer_thickness_update(h,Mak,La_old,input,fid_log,kt)}
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

function La=active_layer_thickness_update(h,Mak,La_old,input,fid_log,kt)
%comment out fot improved performance if the version is clear from github
% version='1';
% if kt==1; fprintf(fid_log,'active_layer_thickness_update version: %s\n',version); end 

%%

switch input.mor.Latype
    case 1 %constant active layer thickness
%         La=La_old; %incorrect if it is changed due to ellipticity solution
        La=repmat(input.mor.La,1,input.mdv.nx);   
    case 2 %active layer thickness related to grain size
        error('not yet implemented')
    case 3 %active layer thickness related to flow depth
        error('not yet implemented')      
    case 4 %active layer thickness growing with time
        La=La_old+input.mor.La_t_growth*input.mdv.dt;
end
        