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
%$Id: ini_Cfb.m 106 2017-06-27 12:45:13Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/ini_Cfb.m $
%
%ini_Cfb is a function that creates the initial bed friction coefficient vector
%
%Cf=ini_Cf(input,fid_log)
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

function Cfb=ini_Cfb(input,fid_log)

%% RENAME

nx=input.mdv.nx;

%% CALC

% 0=NO; 1=Johnson (1942); 2=Nikuradse; 3=imposed bed friction coefficient
switch input.frc.wall_corr
    case 0
        Cfb=input.mdv.Cf(1).*ones(1,nx);              
    case 1 
        error('this should be done...')
    case 2
        error('this should be done...')
    case 3
        Cfb=input.frc.Cfb(1).*ones(1,nx);              
end

end %function
