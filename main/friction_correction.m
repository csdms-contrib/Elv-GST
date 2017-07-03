%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 106 $
%$Date: 2017-06-27 14:45:13 +0200 (Tue, 27 Jun 2017) $
%$Author: V $
%$Id: friction_correction.m 106 2017-06-27 12:45:13Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/friction_correction.m $
%
%friction_correction corrects the friction coefficient for wall friction and ripples
%
%[Cf_b]=friction_correction(u,h,Cf,input,fid_log,kt)
%
%INPUT:
%   -
%
%OUTPUT:
%   -
%
%HISTORY:
%170324
%   -V. Created for the first time.
%


function [Cf_b]=friction_correction(u,h,Cf,Mak,La,input,fid_log,kt)


%% 
%% RENAME
%% 

g=input.mdv.g;
nx=input.mdv.nx;
dk=input.sed.dk;
nf=input.mdv.nf;

%%
%% CALC
%%

%wall correction
switch input.frc.wall_corr %flume wall correction:   ; 
    case 0 %0=NO
        Cf_b=Cf;
        u_st_b=Cf_b.*u.^2;        
    case 1 %1=Johnson (1942)
        Sf=Cf.*u.^2/g./h; %friction slope
        Cf_b=NaN(1,nx); %bed friction coefficient
        u_st_b=NaN(1,nx); %bed friction velocity
        for kx=1:nx
            [~,Cf_b(kx),~,u_st_b(kx)]=side_wall_correction(input,u(kx),h(kx),Sf(kx));
        end
    case 2 %2=Nikuradse
        if nf==1 %unisize calculation (treated differently due to effective fractions)
            Fak=ones(nx,1);
        else %multisize calculation
            Fak=NaN(nx,nf); %preallocate volume fractions
            Fak(:,1:nf-1)=Mak./repmat(La,1,nf-1); %effective fractions
            Fak(:,nf)=ones(nx,1)-sum(Fak(:,1:nf-1),2); %all fractions
        end %isempty(Mak)
        switch input.tra.Dm
            case 1 %geometric
                Dm=2.^sum(Fak.*repmat(log2(dk),nx,1),2);
            case 2 %arithmetic
                Dm=sum(Fak.*repmat(dk,nx,1),2);
        end
        Cf_b=1./(5.75*log10(12*h./(2*Dm))).^2;
    case 3 %3=imposed bed friction coefficient
        Cf_b=repmat(input.frc.Cfb,nx,1); 
    otherwise 
        error('ups... this makes no sense')
end

%ripple correction
if input.frc.ripple_corr==1 %ripple correction: 0=NO; 1=YES
    for kx=1:nx
        [Cf_b(kx),~,~,~]=bed_form_correction(input,u_st_b(kx),u(kx));
    end
end


