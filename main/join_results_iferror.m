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
%$Id: join_results_iferror.m 106 2017-06-27 12:45:13Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/join_results_iferror.m $
%
%joint_results_iferror jpint the temporary results in the output file in case there has been an error.
%
%\texttt{joint_results_iferror(path_file_input,fid_log)}
%
%INPUT:
%   -\texttt{path_file_input} = path to the input.mat file [-]; [(nf-1)x(nx) double]
%   -\texttt{fid_log} = identificator of the log file
%
%OUTPUT:
%   -
%
%HISTORY:
%170531
%   -V. Created for the first time.
%


function joint_results_iferror(path_file_input,fid_log)

%%
%% RENAME
%%

input=NaN; %V is stupid and has just realised that 'input' is also a function in MatLab. GEFELICITEERD!
load(path_file_input); %(input)
        
%%
%% CALL
%%

switch input.mdv.savemethod
    case 2
        join_results(input,fid_log)
end


end %function