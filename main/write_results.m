%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 81 $
%$Date: 2017-04-20 15:51:02 +0200 (Thu, 20 Apr 2017) $
%$Author: V $
%$Id: write_results.m 81 2017-04-20 13:51:02Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/write_results.m $
%
%write_results does this and that
%
%write_results(input,fid_log,kts)
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

function write_results(input,fid_log,kts)
%comment out fot improved performance if the version is clear from github
% version='1';
% if kts==1; fprintf(fid_log,'write_results version: %s\n',version); end 

%% all in one file direclty

if input.mdv.savemethod==1
output_mat=matfile(input.mdv.path_file_output,'writable',true); %matfile io object creation

for ko=1:input.mdv.no
    
    %call to the workspace of the function the variables that we want to save
    aux_varname=input.mdv.output_var{1,ko}; %variable name to update in output.mat
    aux_var=evalin('caller',aux_varname); %variable value in the main function corresponding to the variable name
    feval(@()assignin('caller',aux_varname,aux_var)) %rename such that the variable name goes with its variable value
    
    %update value
    nel=size(output_mat.(aux_varname)); %size of the variable in the .mat file
    if strcmp(aux_varname,'time_loop')
%         time_loop_tmp=eval(aux_varname);
        output_mat.(aux_varname)((kts-2)*input.mdv.Flmap_dt/input.mdv.dt+1:(kts-1)*input.mdv.Flmap_dt/input.mdv.dt,1)=time_loop;
    else
        output_mat.(aux_varname)(1:nel(1),1:nel(2),1:nel(3),kts)=eval(aux_varname); %add the value of the current time
    end
    
    %if problems with dimensions check this:
%     eval(sprintf('%s(%i,%i)',aux_varname,1,2)); %add the value of the current time
    
end

%% separate files
else

for ko=1:input.mdv.no
    %call to the workspace of the function the variables that we want to save
    aux_varname=input.mdv.output_var{1,ko}; %variable name to update in output.mat
    aux_var=evalin('caller',aux_varname); %variable value in the main function corresponding to the variable name
    feval(@()assignin('caller',aux_varname,aux_var)) %rename such that the variable name goes with its variable value
end
	
path_file_output_sng=fullfile(input.mdv.path_folder_TMP_output,sprintf('%06d.mat',kts));
save(path_file_output_sng,input.mdv.output_var{1,:},'-v6') %v6 is faster than v7.3

end
