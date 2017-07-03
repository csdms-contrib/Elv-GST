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
%$Id: run_ELV.m 49 2017-04-04 07:33:12Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0058/main/run_ELV.m $
%
%get_bouss_coef does this and that
%
%alpha_b = get_bouss_coef(input,x,hh)
%
%INPUT:
%   -
%
%OUTPUT:
%   -
%
%HISTORY:
%170404
%   -V. Added the header because Liselot does not follow the protocol :D
%

function alpha_b = get_bouss_coef(input,x,hh)
% Compute the coefficient 
% The computation is performed in different parts for the main channel and
% flood planes seperately.
%
% x should be a node identifier, and h a vector with h(1,x) the location at x!


%% Check wether flow depth is a vector or a single value
if numel(hh)==1
    h = hh;
else
    h = hh(1,x);
end

%% Get cross_sectional_area pars
[Af, Af1, Af2, Af3] = get_cross_section(input,x,h,'Af');
[~, P1, P2, P3] = get_cross_section(input,x,h,'P');


%% Contribution to parameters for each of the subparts
C1 = sqrt(input.mdv.g/input.mdv.Cf(1,1));
C2 = sqrt(input.mdv.g/input.mdv.Cf(2,1));
C3 = sqrt(input.mdv.g/input.mdv.Cf(3,1));
C_v = [C1,C2,C3]; %Chezy values
Af_v = [Af1,Af2,Af3];
P_v = [P1,P2,P3];
R_v = Af_v./P_v;
R_v(isnan(R_v)) = 0;

alpha_b = sum(C_v.^2.*Af_v.*R_v)*Af/(sum(C_v.*Af_v.*sqrt(R_v))^2);

end