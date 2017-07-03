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
%$Id: output_creation.m 106 2017-06-27 12:45:13Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/output_creation.m $
%
%output_creation is a function that creates the file where the results will be saved and saves the initial condition.
%
%output_creation(input,fid_log)
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

function output_creation(input,fid_log)
%comment out fot improved performance if the version is clear from github
% version='1';
% fprintf(fid_log,'output_creation version: %s\n',version); 

%%
%% MAP VARIABLES
%%

%the variables in 'ini_var' are created as initial condition. If more
%variables are asked to be outputed they need to also be created or when
%saving the initial condition and creating the output variable there will
%be an error because the variable is inexisting.
ini_var={'u','h','etab','Mak','La','msk','Ls','Cf'};
for ko=1:input.mdv.no
    if any(strcmp(input.mdv.output_var{1,ko},ini_var)) %if a variabke that you want to save is one of the initial condition ones
        aux_varname=input.mdv.output_var{1,ko}; %variable name to update in output.mat
        aux_var=evalin('caller',aux_varname); %variable value in the main function corresponding to the variable name
        aux_var_size=size(aux_var);
        if numel(aux_var_size)==2
            aux_var_4=NaN(aux_var_size(1),aux_var_size(2),1,input.mdv.nT); 
        else 
            aux_var_4=NaN(aux_var_size(1),aux_var_size(2),aux_var_size(3),input.mdv.nT); 
        end
        aux_var_4(:,:,:,1)=aux_var;
        feval(@()assignin('caller',aux_varname,aux_var_4)) %rename such that the variable name goes with its variable value
    else
        %it is better for plotting purposes to preallocate with the right size if known.
        %qbk
        switch input.mdv.output_var{1,ko}
            case 'qbk'
                qbk=NaN(input.mdv.nf,input.mdv.nx,1,input.mdv.nT); %#ok
            case 'ell_idx'
                ell_idx=false(1,input.mdv.nx,1,input.mdv.nT); %#ok
            case 'time_loop'
                time_loop=NaN(input.mdv.nt,1); %#ok
            case 'pmm'
                pmm=NaN(2,input.mdv.nx,1,input.mdv.nT); %#ok [(alpha,beta),nx,1,nT]
            otherwise
                eval(sprintf('%s=NaN(%d,%d,%d,%d)',input.mdv.output_var{1,ko},nf,nx,nsl,input.mdv.nT));
        end
        
        
    end
end

%% SAVE

%the variables cannot be saved in a structure if you want to partially load
%them. the file must be save in version 7.3 for partially loading it
path_file_output=input.mdv.path_file_output;
save_var_str=input.mdv.output_var; 
save_var_str{input.mdv.no+1}='fieldNames'; %to use v2struct it is necessary to have the string 'fieldNames' in the cell array.
save_var_sctruct=v2struct(save_var_str); %#ok variables into a structure to save it
save(path_file_output,'-struct','save_var_sctruct','-v7.3')

