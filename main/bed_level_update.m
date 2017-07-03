%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 108 $
%$Date: 2017-06-27 16:46:40 +0200 (Tue, 27 Jun 2017) $
%$Author: V $
%$Id: bed_level_update.m 108 2017-06-27 14:46:40Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/bed_level_update.m $
%
%bed_level_update updates the bed elevation
%
%etab_new=bed_level_update(etab,qbk,bc,input,fid_log,kt)
%
%INPUT:
%   -input = variable containing the input [struct] e.g. input
%
%OUTPUT:
%   -
%
%HISTORY:
%160223
%   -V. Created for the first time.
%
%160429
%   -V. Introduction of periodic boundary conditions
%
%160623
%   -V. Cyclic boundary conditions.
%
%160803
%	-L. Merging; including cycled boundary conditions
%
%170126
%   -L. Added cases 13,14 (no new version)
%
%170516
%   -V. Erased upwind factor


function etab_new=bed_level_update(etab,qbk,bc,input,fid_log,kt)

%%
%% RENAME
%%

dx=input.grd.dx;
dt=input.mdv.dt;    
MorFac=input.mor.MorFac;
cb=1-input.mor.porosity;
nx=input.mdv.nx; %number of cells
UpwFac=input.mdv.UpwFac;

%%
%% BOUNDARY CONDITION
%%

switch input.bcm.type
%%upstream bed load
    case {1,3,11,12,13,14,'set1'}
        Qbkt = mod(kt,bc.repQbkT(2))+(mod(kt,bc.repQbkT(2))==0)*bc.repQbkT(2);
        qbk0=input.grd.B(1)*bc.qbk0(Qbkt,:); %[1xnf double]
        qb0=sum(qbk0,2); %[1x1 double]
    case 2
        qbk0=qbk(:,1);
        qb0=sum(qbk0,1); %[1x1 double]       
    otherwise
        error('Kapot! check input.bcm.type')
end

%%
%% EXNER
%%

switch input.mor.bedupdate
    case 0
        etab_new=etab;
        
    case 1
        switch input.mdv.flowtype
            case {1,6}
                %total load
                qb=input.grd.B.*sum(qbk,1); %[1xnf double]
                etab_new(1,1) = etab(1,1) - MorFac * dt./cb.* ((UpwFac * ((qb(1)-qb0)./(dx/2)) + (1-UpwFac) * ((qb(2)-qb(1))./(dx/2)))./input.grd.B(1));
                etab_new(1,nx) = etab(1,nx) - MorFac * dt/cb * ((qb(nx)-qb(nx-1))/(dx))/input.grd.B(end);  
                etab_new(1,2:nx-1) = etab(1,2:nx-1) - MorFac * dt./cb.* ((UpwFac * ((qb(2:nx-1)-qb(1:nx-2))./(dx)) + (1-UpwFac) * ((qb(3:nx)-qb(2:nx-1))./(dx)))./input.grd.B(3:nx));
             
            case {3,4}
                UpwFac = 1-(qbk<0); %sets the UpwFac to 1 if flow comes from left, and to 0 if flow comes from right

                %total load
                qb=input.grd.B.*sum(qbk,1); %[1xnf double]
                etab_new(1,1) = etab(1,1) - MorFac * dt./cb.* ((UpwFac(1) * ((qb(1)-qb0)./(dx/2)) + (1-UpwFac(1)) * ((qb(2)-qb(1))./(dx/2)))./input.grd.B(1));
                
                if qbk(nx)>0              
                    etab_new(1,nx) = etab(1,nx) - MorFac * dt/cb * ((qb(nx)-qb(nx-1))/(dx))/input.grd.B(end);
                else
                    etab_new(1,nx) = etab(1,nx);
                end
                          
                etab_new(1,2:nx-1) = etab(1,2:nx-1) - MorFac * dt./cb.* (UpwFac(2:nx-1).* ((qb(2:nx-1)-qb(1:nx-2))./(dx)) + (1-UpwFac(2:nx-1)).* ((qb(3:nx)-qb(2:nx-1))./(dx)))./input.grd.B(2:end-1);
                
            otherwise
                error('Supposedly you do not end up here');

        end
    otherwise
       error(':( I thought you wanted to use Exner... :(')
        
end
end
