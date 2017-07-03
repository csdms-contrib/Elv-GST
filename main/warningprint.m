%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 107 $
%$Date: 2017-06-27 14:56:45 +0200 (Tue, 27 Jun 2017) $
%$Author: V $
%$Id: warningprint.m 107 2017-06-27 12:56:45Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/warningprint.m $
%
%warningprint does this and that
%
%warningprint(fid_log, text)
%
%INPUT:
%   -input = input structure
%
%OUTPUT:
%   -
%
%HISTORY:

function warningprint(fid_log, text)
if isnan(fid_log)==1
    warning('No log file specified, warning is only printed to the screen');
    warning(text);
else
    %   Display warning message
    %   Write output to logfile
    warning(text);
    fprintf(fid_log,'!!WARNING!!: %s \n',text);
end
end

