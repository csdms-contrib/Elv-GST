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
%dv_cross_aux does this and that
%
%diff = dv_cross_aux(input, x, term)
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

function diff = dv_cross_aux(input, x, term)

switch term
    case 'B_main'
        var = input.grd.B(1,:);
    case 'B_l'
        var = input.grd.Bparam(2,:);
    case 'B_r'
        var = input.grd.Bparam(4,:);
    case 'Bfl'
        var = input.grd.B(2,:);
    case 'Bfr'
        var = input.grd.B(3,:);
    case 'h_l'
        var = input.grd.Bparam(1,:);
    case 'h_r'
        var = input.grd.Bparam(3,:);
end

%Add dummy cells;
%var = [var(1)-(var(2)-var(1)), var, var(end)+(var(end)-var(end-1))];

% Compute derivative
%try %centered
if x==1
    diff= (var(x+1)-var(x))/(input.grd.dx);
elseif x==numel(var)
    diff= (var(x)-var(x-1))/(input.grd.dx);
else
    diff= (var(x+1)-var(x-1))/(2*input.grd.dx);
end
%catch 
%    try %upstream 
%        diff= (var(x+1)-var(x))/(input.grd.dx);
%    catch
%        try %downstream
%            diff= (var(x)-var(x-1))/(input.grd.dx);
%        catch
%            %it should be a fixed with case
%            error('check dimensions of width');
%        end
%    end
%end      
end