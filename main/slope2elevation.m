%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 14 $
%$Date: 2017-02-06 16:28:08 +0100 (Mon, 06 Feb 2017) $
%$Author: V $
%$Id: slope2elevation.m 14 2017-02-06 15:28:08Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/slope2elevation.m $
%
%slope2elevation computes the bed elevation given the slope

%INPUT:
%   -input = variable containing the input [struct] e.g. input
%
%OUTPUT:
%   -
%
%HISTORY:
%160223
%   -V. Created for the first time.

function etab=slope2elevation(slopeb,etab0,input,fid_log)
%comment out fot improved performance if the version is clear from github
% version='1';
% fprintf(fid_log,'slope2elevation version: %s\n',version);

%%

etab=NaN(1,input.mdv.nx);
etab(end)=etab0+slopeb(end)*input.grd.dx/2;
for kx=input.mdv.nx-1:-1:1
    etab(kx)=etab(kx+1)+slopeb(kx)*input.grd.dx;
end