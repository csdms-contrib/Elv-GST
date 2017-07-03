%run in the folder where the script is

%% PREAMBLE

fclose all;
clear;
close all;
clc;

%% INPUT 

runid_serie='trial';
runid_number='003';
input_filename='input_ELV';
paths_runs='c:\Users\victorchavarri\temporal\ELV\';
erase_previous=0; %it is dangerous, use with care and attention
do_profile=0; %0=NO; 1=YES
do_postprocessing=0; %0=NO; 1=YES


%%  %%  %%
  %%  %%       DO NOT TOUCH FROM HERE ON
%%  %%  %%    


%% PATH DEFINITIONS

source_path = pwd; 
paths_main              = fullfile(source_path,'..',filesep,'main');
paths_auxiliary         = fullfile(source_path,'..',filesep,'auxiliary');
paths_postprocessing    = fullfile(source_path,'..',filesep,'postprocessing');
paths_input_script      = fullfile(source_path,input_filename);
paths_run_folder        = fullfile(paths_runs,runid_serie,runid_number);
paths_input_mat         = fullfile(paths_run_folder,'input.mat');

%% ERASE EVERYTHING IN THE FOLDER

if erase_previous
    time_d=5; %#ok seconds delay
    warning('off','backtrace');
    for kt=time_d:-1:0
        warning('You are going to erase the simulation %s%s in %d seconds',runid_serie,runid_number,kt);
        pause(1)
    end
    warning('on','backtrace');
    fclose all;
    dire_in=paths_run_folder;
    dire=dir(dire_in);
    for kf=3:numel(dire)
        if ispc
            if exist(fullfile(dire_in,dire(kf).name),'dir')==7
                dos(sprintf('RD /S /Q %s',fullfile(dire_in,dire(kf).name)));
            elseif exist(fullfile(dire_in,dire(kf).name),'file')
                dos(sprintf('DEL %s',fullfile(dire_in,dire(kf).name)));
            end
        elseif ismac
            error('Are you seriously using a Mac? Come on... You will have to manually erase the folder and set erase_previous to 0')
        else
            error('You will have to manually erase the folder and set erase_previous to 0')
        end
    end
end

%% CREATE INPUT

run(paths_input_script);
input.run=strcat(runid_serie,runid_number); %simulation name [char]; e.g. 'L\_05'

%check before saving
if exist(paths_run_folder,'dir')==7
    if exist(paths_input_mat,'file')
        error('It seems you are trying to overwrite results. It already exist a file input.mat in the target folder')
    else
        save(paths_input_mat,'input')
    end
else
    mkdir(paths_run_folder)
    save(paths_input_mat,'input')
end

%% ADD PATHS

%paths to add if they are not already added (needs to be here because the
%run folder is created in the previous step)
paths2add{1,1}=paths_main;
paths2add{2,1}=paths_run_folder;       
paths2add{3,1}=paths_auxiliary;
paths2add{4,1}=paths_postprocessing;

paths_inmatlab=regexp(path,pathsep,'split');
for kp=1:numel(paths2add)
    if ispc  % Windows is not case-sensitive
      onPath=any(strcmpi(paths2add{kp,1},paths_inmatlab));
    else
      onPath=any(strcmp(paths2add{kp,1},paths_inmatlab));
    end
    if onPath==0
        addpath(paths2add{kp,1});
    end
end

%% RUN 

switch do_profile
    case 1
        paths_profile=fullfile(paths_run_folder,'profile');
        mkdir(paths_profile)
        profile on
        run_ELV(paths_input_mat);
        profile off
        profsave(profile('info'),paths_profile)
        profile viewer
    otherwise
        run_ELV(paths_input_mat);
end

%% POST
if do_postprocessing
    %% patch
input_fig_input; %#ok
fig_patch(paths_run_folder,fig_input)
    %% level
input_fig_input;
fig_level(paths_run_folder,fig_input)
end