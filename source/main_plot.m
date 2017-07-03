
%%
clear all %if changes outside matlab
%%
cd d:\victorchavarri\SURFdrive\projects\00_codes\ELV\branch_V\source\

%% INPUT

% dire_in='c:\Users\victorchavarri\temporal\ELV\K\041\'; 
% dire_in='c:\Users\victorchavarri\temporal\ELV\J\042\'; 
dire_in='c:\Users\victorchavarri\temporal\ELV\trial\001\';
% dire_in='d:\victorchavarri\SURFdrive\projects\02_runs\ELV\J\002\';

%% RUN

run('input_fig_input.m')
addpath('..\postprocessing\')

%% patch
fig_patch(dire_in,fig_input)    
%% level
fig_level(dire_in,fig_input)    
%% x-cnt
fig_x_cnt(dire_in,fig_input)    
%% x-t
fig_xt(dire_in,fig_input)    
%% time_loop
fig_time_loop(dire_in,fig_input)