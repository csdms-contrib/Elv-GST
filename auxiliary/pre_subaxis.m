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
%$Id: pre_subaxis.m 7 2017-02-06 13:42:31Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/auxiliary/pre_subaxis.m $
%
% npr=1; %number of plot rows
% npc=1; %number of plot columns
% marg.mt=1.0; %top margin [cm]
% marg.mb=2.0; %bottom margin [cm]
% marg.mr=4.0; %right margin [cm]
% marg.ml=2.0; %left margin [cm]
% marg.sh=0.0; %horizontal spacing [cm]
% marg.sv=0.0; %vertical spacing [cm]
% 
% [mt,mb,mr,ml,sh,sv]=pre_subaxis(h.fig,marg.mt,marg.mb,marg.mr,marg.ml,marg.sh,marg.sv);
% h.sfig=subaxis(npr,npc,1,1,1,1,'mt',mt,'mb',mb,'mr',mr,'ml',ml,'sv',sv,'sh',sh);



function [mt,mb,mr,ml,sh,sv]=pre_subaxis(handle,mt_i,mb_i,mr_i,ml_i,sh_i,sv_i)

dim=get(handle,'paperposition');
width=dim(3);
height=dim(4);

mt=mt_i/height;
mb=mb_i/height;
mr=mr_i/width;
ml=ml_i/width;
sh=sh_i/width;
sv=sv_i/height;

