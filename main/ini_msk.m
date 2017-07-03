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
%$Id: ini_msk.m 79 2017-04-20 11:07:54Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/ini_msk.m $
%
%ini_msk is a function that creates the substrate variable
%
%[msk,Ls]=ini_msk(La,input,fid_log)
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
%160404
%   -V. Patch addition.
%
%160415
%   -V. Bug, Ls had size [nef,nx,nsl] but it should be [1,nx,nsl]. No
%   problems or 1 or 2 fractions. 
%
%170317
%   -V. Add last thick layer
%

function [msk,Ls]=ini_msk(La,input,fid_log)
%comment out fot improved performance if the version is clear from github
% version='3';
% fprintf(fid_log,'ini_msk version: %s\n',version);

%%
%% RENAME
%%

ThUnLyr=input.mor.ThUnLyr;
nef=input.mdv.nef;
nx=input.mdv.nx;
fsk=input.ini.fsk;
nsl=input.mdv.nsl;
ThUnLyrEnd=input.mor.ThUnLyrEnd;

if isfield(input.ini,'subs') 
    if isfield(input.ini.subs,'patch') 
        xedg=input.mdv.xedg;
        patch_x=input.ini.subs.patch.x; 
        patch_releta=input.ini.subs.patch.releta;
        patch_fsk=input.ini.subs.patch.fsk; 
    end
end

%%
%% GENERAL
%%

if size(fsk,2)==1 %if no variation in streamwise direction
    fsk=repmat(fsk,1,nx,nsl);
else 
    fsk=repmat(fsk,1,1,nsl);
end
Ls=ThUnLyr.*ones(1,nx,nsl);
Ls(:,:,end)=ThUnLyrEnd;
Ls(:,:,1)=Ls(:,:,1)-1e-8; %if it is exactly the same length as ThUnLyr it gives problems when calculating cell jumps for the first time. then it is never equal to ThUnLyr.


%% 
%% PATCH
%%

if isfield(input.ini,'subs') 
    if isfield(input.ini.subs,'patch') 

        %x coordinate
        [~,xc(1)]=min(abs(xedg-patch_x(1)));
        [~,xc(2)]=min(abs(xedg-patch_x(2)));
        nxcp=xc(2)-xc(1)+1; %number of xcells in the patch

        %y coordinate
        Lt=NaN(1,nx,nsl+1); %total thickness
        Lt(1,:,1)=La;
        Lt(1,:,2:end)=Ls; 
        Lcum=cumsum(Lt,3);
        [~,yc]=min(abs(patch_releta-Lcum(1,xc(1),:)));
        nycp=nsl-yc+1;

        %modify substrate
        fsk(:,xc(1):xc(2),yc:end)=repmat(patch_fsk,1,nxcp,nycp);    

    end
end

%volume fraction to mass
msk=fsk.*repmat(Ls,nef,1,1);