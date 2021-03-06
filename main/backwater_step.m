%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 14 $
%$Date: 2017-02-06 16:28:08 +0100 (ma, 06 feb 2017) $
%$Author: V $
%$Id: flow_update.m 14 2017-02-06 15:28:08Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/L2/main/flow_update.m $
%
%backwater_step computes just one step of the solver u,h,etab,Cf should be single values;
%
%[U,H]=backwater_step(x,h,ib,Cf,Q,input)
%
%INPUT:
%   -
%
%OUTPUT:
%   -
%
%HISTORY:

function [U,H] = backwater_step(x,h,ib,Cf,Q,input)

[Af,Af1,Af2,Af3] = get_cross_section(input,x,h,'Af');
[P,P1,P2,P3] = get_cross_section(input,x,h,'P');
[b1,b2,~,~,~,~] = dv_cross_section(input,x,h,'Af');
alpha_b = get_bouss_coef(input,x,h);
[a1,a2] = dv_bouss_coef(input,x,h);
R = Af1/P1+Af2/P2+Af3/P3;%Af./P;

%R = h; 
%a1 = 0;
%a2 = 0;
%b2 = 0;
H = h - input.grd.dx*(b2*alpha_b*Q^2/Af^2-a2*Q^2/Af+input.mdv.g*Af*ib-Cf*abs(Q)*Q/(R*Af))/(-b1*alpha_b*Q^2/Af^2+a1*Q^2/Af + input.mdv.g*Af);
%H = h - input.grd.dx*(ib-Cf*abs(Q)*Q/(input.mdv.g*h^3*input.grd.B(1,1)^2))/(1-Q^2/input.grd.B(1,1)^2/h^3);
U = get_flow_velocities(input,x,Q,H,ib);

end
