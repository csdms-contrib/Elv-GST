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
%$Id: active_layer_mass_update.m 108 2017-06-27 14:46:40Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/active_layer_mass_update.m $
%
%active_layer_mass_update updates the mass (volume) of sediment at the active layer.
%
%\texttt{Mak_new=active_layer_mass_update(Mak,detaLa,fIk,qbk,bc,input,fid_log,kt)}
%
%INPUT:
%   -\texttt{Mak} = effective volume of sediment per unit of bed area in the active layer [m]; [(nf-1)x(nx) double]
%   -\texttt{detaLa} = variation in elevation of the interface between the active layer and the substarte [m]; [(1)x(nx) double]
%   -\texttt{fIk} = effective volume fraction content of sediment at the interface between the active layer and the substrate [-]; [(nf-1)x(nx) double]
%   -\texttt{qbk} = volume of sediment transported excluding pores per unit time and width and per size fraction [m^2/s]; [(nf)x(nx) double]
%   -\texttt{bc} = boundary conditions structure 
%   -\texttt{input} = input structure
%   -\texttt{fid_log} = identificator of the log file
%   -\texttt{kt} = time step counter [-]; [(1)x(1) double]
%
%OUTPUT:
%   -\texttt{Mak_new} = new effective volume of sediment per unit of bed area in the active layer [m]; [(nf-1)x(nx) double]
%
%HISTORY:
%160223
%   -V. Created for the first time.
%
%160429
%   -V. Periodic boundary conditions.
%
%160623
%   -V. Cyclic boundary conditions.
%
%160825
%   -L. Repeated hydrograph - no version update
%
%170517
%   -V. Removed upwind factor.

function Mak_new=active_layer_mass_update(Mak,detaLa,fIk,qbk,bc,input,fid_log,kt)

%%
%% RENAME
%%

dx=input.grd.dx;
dt=input.mdv.dt;    
MorFac=input.mor.MorFac;
nef=input.mdv.nef;
cb=1-input.mor.porosity;
nx=input.mdv.nx; 
UpwFac=input.mdv.UpwFac;

%%
%% BOUNDARY CONDITION
%%

switch input.bcm.type
%%upstream bed load
    case {1,3,11,12,13,'set1'}
        Qbkt = mod(kt,bc.repQbkT(2))+(mod(kt,bc.repQbkT(2))==0)*bc.repQbkT(2);
        qbk0=input.grd.B(1)*bc.qbk0(Qbkt,:)'; %[1xnf double]
    case 2
        qbk0=qbk(:,end);  %[1xnf double]  
    otherwise
        error('input.bcm.type...')
end

%%
%% ACTIVE LAYER EQUATION
%%

qbk=repmat(input.grd.B,nef+1,1).*qbk; %[1xnf double]
detaLa_dt=detaLa/dt; %%d(\eta-La)/dt [1xnx double]

Mak_new(:,1)=Mak(:,1)-dt*(fIk(:,1).*detaLa_dt(1,1)+MorFac/cb*((qbk(1:nef,1)-qbk0(1:nef,1))/(dx/2))./input.grd.B(2));
Mak_new(:,2:nx-1)=Mak(:,2:nx-1)-dt*(fIk(:,2:nx-1).*repmat(detaLa_dt(1,2:nx-1),nef,1)+MorFac/cb*(UpwFac*(qbk(1:nef,2:nx-1)-qbk (1:nef,1:nx-2))/dx+(1-UpwFac)*(qbk(1:nef,3:nx)-qbk(1:nef,2:nx-1))/dx)./(repmat(input.grd.B(3:nx),nef,1)));%input.grd.B(3:K)');
Mak_new(:,nx    )=Mak(:,nx    )-dt*(fIk(:,nx    ).*       detaLa_dt(1,nx    )       +MorFac/cb*(       (qbk(1:nef,nx    )-qbk (1:nef,nx-1  ))/dx)./input.grd.B(nx));


