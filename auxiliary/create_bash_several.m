%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 49 $
%$Date: 2017-04-04 09:33:12 +0200 (Tue, 04 Apr 2017) $
%$Author: V $
%$Id: create_bash_several.m 49 2017-04-04 07:33:12Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/auxiliary/create_bash_several.m $
%
%this function creates a bash file to be run by the cluster

%INPUT:
%   -
%
%OUTPUT:
%   -
%
%HISTORY:
%170403
%   -V. Created for the first time.

function create_bash_several(paths_bash,paths_serie,paths_source)
%% RENAME


%% FILE

data{1  ,1}=        '#!/bin/bash';
data{2  ,1}=sprintf('maindir=%s',paths_serie);
data{3  ,1}=sprintf('exedir=%s' ,fullfile(paths_source,'SINGLE_RUN'));
data{4  ,1}=sprintf('exedirm=%s',fullfile(paths_source,'single_run_ELV.m'));
data{5  ,1}=        'for dir in $maindir/*/';
data{6  ,1}=        'do';
data{7  ,1}=        '	cd $dir';
data{8  ,1}=        '	cp $exedir SINGLE_RUN';
data{9  ,1}=        '	cp $exedirm single_run_ELV.m';
data{10 ,1}=        '';
data{11 ,1}=        '	job_name=PARS_$dir';
data{12 ,1}=        '	qsub -N $job_name  SINGLE_RUN';
data{13 ,1}=        'done';

%% WRITE

file_name=paths_bash;

%check if the file already exists
if exist(file_name,'file')
    error('You are trying to overwrite a file!')
end

fileID_out=fopen(file_name,'w');
for kl=1:numel(data)
    fprintf(fileID_out,'%s\n',data{kl,1}); %attention, no \r or unix will complain
end

fclose(fileID_out);