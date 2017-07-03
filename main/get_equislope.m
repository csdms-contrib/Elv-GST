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
%$Id: get_equislope.m 79 2017-04-20 11:07:54Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/get_equislope.m $
%
%get_equislope does this and that
%
%S = get_equislope(Qw,input,AL,Fk,X0)
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
%

function S = get_equislope(Qw,input,AL,Fk,X0)
    if nargin<5
        X0 = 1e-4;
    end
    %input = add_sedflags(input);
    F = @(X)solve_nfbc([X Fk(1:end-1)],input,Qw,AL);
    S = fzero(F,X0);
    disp('Objective value of total sediment load (should be zero):');
    solve_nfbc([S Fk(1:end-1)],input,Qw,AL)
end
