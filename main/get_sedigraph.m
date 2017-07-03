%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 79 $
%$Date: 2017-04-20 13:07:54 +0200 (Thu, 20 Apr 2017) $
%$Author: V $
%$Id: get_sedigraph.m 79 2017-04-20 11:07:54Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/get_sedigraph.m $
%
%get_sedigraph does this and that
%
%Qb = get_sedigraph(Qw,input,S,Fk)
%
%INPUT:
%   -input = input structure
%
%OUTPUT:
%   -
%
%HISTORY:
%161128
%   -L. Created for the first time

function Qb = get_sedigraph(Qw,input,S,Fk)
% S is equilibrium slope
% Fk is all the fractions; we select all but the coarsest which is
% 1-sum(others)
input = add_sedflags(input);
if numel(input.sed.dk)>1
    La = ones(1,numel(Qw));
    Mak = repmat(Fk(1:end-1)',1,numel(Qw));
else
    La = NaN*ones(1,numel(Qw));
    Mak = NaN*ones(1,numel(Qw));
end    
h = (input.mdv.Cf(1).*(Qw/input.grd.B(1,1)).^2/(9.81*S)).^(1/3);
Cf = input.mdv.Cf(1).*ones(1,numel(Qw)); 
if isfield(input,'tra.calib')==1
else 
    input.tra.calib =1;
end
[qbk,~]=sediment_transport(input.aux.flg,input.aux.cnt,h,Qw/input.grd.B(1,1),Cf,La,Mak',input.sed.dk,input.tra.param,input.aux.flg.hiding_parameter,input.tra.calib,NaN,NaN)
Qb = qbk*input.grd.B(1,1);
end
