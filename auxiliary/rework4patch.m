%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 7 $
%$Date: 2017-02-06 14:42:31 +0100 (Mon, 06 Feb 2017) $
%$Author: V $
%$Id: rework4patch.m 7 2017-02-06 13:42:31Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/auxiliary/rework4patch.m $
%
function [f,v,col]=rework4patch_160302(in)

%% RENAME in

XCOR=in.XCOR;

sub=in.sub;
cvar=in.cvar;

nx=in.nx;
nl=in.nl;

%%

ncx=nx-2;
ncordx=nx-1;
ncordy=nl+1;

%% data rework

%patches
Cm=cvar'; %fraction 2
Xm=repmat(XCOR,ncordy,1);
Ym=sub;

xr=reshape(Xm',ncordx*ncordy,1);
yr=reshape(Ym',ncx*(nl+1),1);
aux= reshape(repmat(1:1:ncordx,2,1),ncordx*2,1);
aux1=reshape(repmat(1:1:ncx*(nl+1),2,1),(ncx*(nl+1))*2,1);
aux2=aux(2:end-1);
aux3=repmat(aux2,nl+1,1);
v(:,1)=xr(aux3);
v(:,2)=yr(aux1);

f=NaN(nl*ncx,4);
f(:,1:2)=(reshape(1:1:(ncx*2*nl),2,nl*ncx))';
aux=(reshape((ncx*2+1):1:((ncx*2)*(nl+1)),2,nl*ncx))';
f(:,3)=aux(:,2);
f(:,4)=aux(:,1);

col=reshape(Cm,ncx*nl,1);