% This is a function which is only a temporary solution to an annoying
% problem. 

function input = add_sedflags(input)

%% MDV

input.mdv.nf=numel(input.sed.dk); %number of sediment size fraction [-]
input.mdv.nef=input.mdv.nf-1;

%% TRA

flg.sed_trans=input.tra.cr;
flg.friction_closure=1;


if input.mdv.nf~=1 %mixed-sediment
    if isfield(input.tra,'hid') %hiding correction
        flg.hiding=input.tra.hid;    
        if input.tra.hid==2 %Power Law
            flg.hiding_parameter=input.tra.hiding_b; %trick to avoid another variable in main function
        else 
            flg.hiding_parameter=NaN;
        end
    else %no hiding correction
        flg.hiding=0;
        flg.hiding_parameter=NaN; %trick to avoid another variable in main function
    end
    if isfield(input.tra,'param')==0
        input.tra.param=NaN;
    end
    if isfield(input.tra,'Dm')==0 %hiding correction
        flg.Dm=1; %default is geometric
    else
        flg.Dm=input.tra.Dm;
    end
else %unisize
    flg.hiding=0;
    flg.Dm=1; 
    flg.hiding_parameter=NaN; %trick to avoid another variable in main function
end

cnt.g=input.mdv.g;
cnt.rho_w=input.mdv.rhow;
cnt.p=0; %this is to compute the sediment transport without pores. Porosity is in bed level update.
cnt.R=(input.sed.rhos-input.mdv.rhow)/input.mdv.rhow;

input.aux.flg=flg;
input.aux.cnt=cnt;


%input.mdv.nf=numel(input.sed.dk); %number of sediment size fraction [-]
%%% flags and constant for the sediment transport function
%%the sediment transport function requires these structures. This can be improved by changing sediment_transport 
%flg.sed_trans=input.tra.cr;
%flg.friction_closure=1;
%
%if input.mdv.nf~=1 %mixed-sediment
%    if isfield(input.tra,'hid') %hiding correction
%        flg.hiding=input.tra.hid;
%        flg.Dm=1;
%    
%        if input.tra.hid==2 %Power Law
%            flg.hiding_parameter=input.tra.hiding_b; %trick to avoid another variable in main function
%        else 
%            flg.hiding_parameter=NaN;
%        end
%    else %no hiding correction
%        flg.hiding=0;
%        if input.sed.tra==4 %WC requires Dm 
%            flg.Dm=1; 
%        else
%            flg.Dm=NaN;
%        end
%        flg.hiding_parameter=NaN; %trick to avoid another variable in main function
%    end
%    if isfield(input.tra,'param')==0
%        input.tra.param=NaN;
%    end
%else %unisize
%    flg.hiding=0;
%    flg.Dm=NaN; 
%    flg.hiding_parameter=NaN; %trick to avoid another variable in main function
%end
%
%cnt.g=input.mdv.g;
%cnt.rho_w=input.mdv.rhow;
%cnt.p=0;
%cnt.R=(input.sed.rhos-input.mdv.rhow)/input.mdv.rhow;
%
%input.aux.flg=flg;
%input.aux.cnt=cnt;
%
%%UNTITLED Summary of this function goes here
%%   Detailed explanation goes here
%
%
end

