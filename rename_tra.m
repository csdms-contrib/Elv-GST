
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
