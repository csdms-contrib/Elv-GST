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
%$Id: condition_construction.m 107 2017-06-27 12:56:45Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/condition_construction.m $
%
%condition_construction does this and that
%
%[u,h,etab,Mak,La,msk,Ls,Cf,bc] = condition_construction (input,fid_log)
%
%INPUT:
%   -input = input structure
%   -fid_log = identificator of the log file
%
%OUTPUT:
%   -
%
%HISTORY:

function [u,h,etab,Mak,La,msk,Ls,Cf,bc] = condition_construction (input,fid_log)
% try
%     version='2';
%     fprintf(fid_log,'condition_construction: %s\n',version);
% catch
% end

%condition_construction is a function that regulates the order in which the initial and boundary conditions are constructed

%INPUT:
%   -input = variable containing the input [struct] e.g. input
%
%OUTPUT:
%   -bc = boundary conditions [struct] 
%
%HISTORY:
%160429
%   -L. Created for the first time.
%170127
%   -L. Updated case numbers

switch input.ini.initype
    case {1,2,3,4}
        switch input.bch.uptype
            case {1,11}
                % Initial condition
                fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of initial condition construction');
                [u,h,etab,Mak,La,msk,Ls,Cf]=initial_condition_construction(input,fid_log);

                % Boundary condition
                fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of boundary conditions construction');
                bc=boundary_conditions_construction(u,h,Mak,La,Cf,input,fid_log);
            case {12,14}
                % Boundary condition
                fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of boundary conditions construction');
                bc=boundary_conditions_construction(NaN,NaN,NaN,NaN,NaN,input,fid_log);               
                
                input.bcm.Qbk0(1,:) = bc.qbk0(1,:)*input.grd.B(1);
                input.bch.Q0(1) = bc.q0(1)*input.grd.B(1);
                
                % Initial condition
                fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of initial condition construction');
                [u,h,etab,Mak,La,msk,Ls,Cf]=initial_condition_construction(input,fid_log);

                % Boundary condition
                fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of boundary conditions construction');
                bc=boundary_conditions_construction(u,h,Mak,La,Cf,input,fid_log);
        end
        
    case {5,12,13,51,52,53}
        % Boundary condition
        fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of boundary conditions construction');
        bc=boundary_conditions_construction(NaN,NaN,NaN,NaN,NaN,input,fid_log);      
        
        % Initial condition
        fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of initial condition construction');
        [u,h,etab,Mak,La,msk,Ls,Cf]=initial_condition_construction(input,fid_log,bc);
        
end


