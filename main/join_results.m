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
%$Id: join_results.m 106 2017-06-27 12:45:13Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/join_results.m $
%
%join_results does this and that
%
%join_results(input,fid_log)
%
%INPUT:
%   -
%
%OUTPUT:
%   -
%
%HISTORY:
%170207
%   -V. Created for the first time.

function join_results(input,fid_log)

%%
%% RENAME
%%

savemethod=input.mdv.savemethod;

%%

if savemethod==2
fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of putting all result files together');

%load the empty results
output_m=load(input.mdv.path_file_output);

%check how many files are in the temporary output folder (do not use nT in
%case the simulation has crashed)
dir_TMP_output=dir(input.mdv.path_folder_TMP_output);
nF=numel(dir_TMP_output)-2; %number of files in directory (. and ..)
nTt=nF+1; %number of saving results times. 't' because it may has crashed so it is 'temporary' :) we add 1 because file 000001 does not exist.

%load the separate resutls files and copy to the variable with all the results
for kT=2:nTt %loop on separate files with results (kT is counter for results time file)
    path_file_output_sng=fullfile(input.mdv.path_folder_TMP_output,sprintf('%06d.mat',kT)); %path to the separate file with results
    output_par=load(path_file_output_sng); %load the partial results
    for ko=1:input.mdv.no %loop on varaibles to save
        aux_varname=input.mdv.output_var{1,ko}; %variable name to update in output.mat
        switch aux_varname 
            case 'time_loop' %history variable
                output_m.(aux_varname)((kT-2)*floor(input.mdv.Flmap_dt/input.mdv.dt)+1:(kT-1)*floor(input.mdv.Flmap_dt/input.mdv.dt))=output_par.(aux_varname);
            case 'pmm'
                output_m.(aux_varname)(:,:,:,kT)=[output_par.(aux_varname).alpha;output_par.(aux_varname).beta];
            otherwise %map variable
                nel=size(output_m.(aux_varname)); %size of the variable in the .mat file
                output_m.(aux_varname)(1:nel(1),1:nel(2),1:nel(3),kT)=output_par.(aux_varname)(1:nel(1),1:nel(2),1:nel(3));
        end
    end
end
        
%modify the variable with all the results
output_mat=matfile(input.mdv.path_file_output,'writable',true); %matfile io object creation
for ko=1:input.mdv.no
    aux_varname=input.mdv.output_var{1,ko}; %variable name to update in output.mat
    output_mat.(aux_varname)=output_m.(aux_varname); %add the value of the current time
end

%% erase single files

%erase first files and then folder as a safety mechanism
%     dire=dir(input.mdv.path_folder_TMP_output);
%     for kf=3:numel(dire)
%         [~,~,ext]=fileparts(dire(kf).name);
%         if strcmp(ext,'.mat') %safety mechanism, only .mat files
%             dos(sprintf('DEL %s',fullfile(dire(kf).folder,dire(kf).name)));
%         end
%     end

%erase folder directly
if ispc
    dos(sprintf('RD /S /Q %s',input.mdv.path_folder_TMP_output));
elseif isunix
    system(sprintf('rm -rf %s',input.mdv.path_folder_TMP_output));
elseif ismac
    warningprint(fid_log,'Are you seriously using a mac? come on... :( very disappointing...');
else
    warningprint(fid_log,'What kind of operating system are you using? Whatever, I cannot erase the output files! :(');
end

fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'End of putting all result files together');
end
