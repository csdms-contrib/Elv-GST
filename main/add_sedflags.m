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
%$Id: add_sedflags.m 107 2017-06-27 12:56:45Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/add_sedflags.m $
%
%add_sedflags parse the sediment input and adds the sediment transport relation flags to the input structure 
%
%\texttt{input=add_sedflags(input)}
%
%INPUT:
%   -input = input structure
%
%OUTPUT:
%   -input = input structure
%
%HISTORY:
%160223
%   -V. Created for the first time.

function input = add_sedflags(input,fid_log)

if nargin<2
    fid_log =NaN;
end

if isfield(input.tra,'cr')==0 
    input.tra.cr=1;
    warningprint(fid_log,'You have not specified a sediment transport relation (input.tra.cr). I will use Meyer-Peter Muller (1948) just because I like it.')
end

if isfield(input.tra,'param')==0 
    switch input.tra.cr
        case 1
            warningprint(fid_log,'You have not specified the parameters of the sediment transport relation (input.tra.param). I will use the standard values for Meyer-Peter and Muller. Attention! The standard is without hiding.')
            input.tra.param=[8,1.5,0.047];
            input.tra.hid=0;
        case 2 
            warningprint(fid_log,'You have not specified the parameters of the sediment transport relation (input.tra.param). I will use the standard values, good luck!')
            input.tra.param=[0.05,5];
        case 3
            warningprint(fid_log,'You have not specified the parameters of the sediment transport relation (input.tra.param). I will use the standard values for Ashida-Michiue. Attention! The standard is with hiding.')
            input.tra.param=[17,0.05];
            input.tra.hid=3;
    end
end

if isfield(input.tra,'hid')==0 %hiding correction  
    input.tra.hid=0;
end

if isfield(input.tra,'hiding_b')==0 %hiding parameter of power law
    input.tra.hiding_b=NaN;
end

if isfield(input.tra,'Dm')==0
    input.tra.Dm=1; %default is geometric
end

%the sediment transport function requires these structures. This can be improved by changing sediment_transport 
flg.sed_trans=input.tra.cr;
flg.friction_closure=1;
flg.hiding=input.tra.hid;
flg.Dm=input.tra.Dm; 
flg.hiding_parameter=input.tra.hiding_b;

cnt.g=input.mdv.g;
cnt.rho_w=input.mdv.rhow;
cnt.p=0; %this is to compute the sediment transport without pores. Porosity is in bed level update.
cnt.R=(input.sed.rhos-input.mdv.rhow)/input.mdv.rhow;

input.aux.flg=flg;
input.aux.cnt=cnt;

end

