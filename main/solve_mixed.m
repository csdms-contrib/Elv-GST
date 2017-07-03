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
%$Id: solve_mixed.m 107 2017-06-27 12:56:45Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/solve_mixed.m $
%
function obj = solve_mixed(X,input,Fak_old, Qbk, dQbkdu, h_old, Fr_old, u_old, pq, K, AL, dxi, dq)
%solve_mixed computes an update of the space marching algorithm
% VERSION 3

%INPUT:
%   -input
%
%OUTPUT:
%   -

%HISTORY:
%

ib = X(1);
dFak = X(2:end);

c_f = input.mdv.Cf;
nf = input.mdv.nf; 

% flow variables
htemp = h_old - (ib - c_f.*Fr_old.^2)./(1-Fr_old.^2)*dxi;
     
% check total load
Fak_new = Fak_old + dFak*dxi;
Mak_new = repmat(Fak_new.*input.mor.La,K,1);
Cf = input.mdv.Cf.*ones(1,K);
La = ones(1,K);
qbk_new = sediment_transport(input.aux.flg,input.aux.cnt,htemp,dq,Cf,La,Mak_new',input.sed.dk,input.tra.param,input.aux.flg.hiding_parameter,1,NaN,NaN);                    
CL = sum(repmat(pq,1,nf).*qbk_new);

% objective
obj = AL - CL;
end

