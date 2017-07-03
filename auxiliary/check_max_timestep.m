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
%$Id: check_max_timestep.m 49 2017-04-04 07:33:12Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/auxiliary/check_max_timestep.m $
%
% This is a very expensive function to check the maximum time step

function [max_dt_qs,max_dt_fc]=check_max_timestep(input,u,h,Cf,La,Fak,Fik)
%% NO INPUT (necessary for CFL)

input=add_sedflags(input);

flg.sed_trans=input.tra.cr;
flg.read=1;
flg.check_input=1;
flg.friction_input=1;
flg.derivatives=1;
flg.cp=0;
flg.anl=[1,2];
    % 1 = fully coupled
    % 2 = quasi-steady
flg.friction_closure=1;
flg.hiding=input.tra.hid;
flg.Dm=input.tra.Dm;
hiding=input.tra.hiding_b;
sedTrans=input.tra.param;
gsd=input.sed.dk';

v=0;
u_b=0.001;  
I=0.1;
twoDnr=30;

%% ECT

if any(flg.anl==1) || any(flg.anl==2) || any(flg.anl==3) || any(flg.anl==4) || any(flg.anl==5)
nodeState=[u,h,Cf,La,Fak,1-sum(Fak),Fik,1-sum(Fik)]; %1D    
elseif (any(flg.anl==6) || any(flg.anl==7) || any(flg.anl==8) || any(flg.anl==9) || any(flg.anl==10) || any(flg.anl==11))
nodeState=[u,v,h,Cf,La,Fak,1-sum(Fak),Fi1,1-sum(Fik)]; %2D
end

[eigen_all,elliptic,A,cp,Ribb,out,...
  eigen_all_qs,elliptic_qs,A_qs,...
  eigen_all_dLa,elliptic_dLa,A_dLa,...
  eigen_all_ad,elliptic_ad,A_ad,...
  eigen_all_2Dx,eigen_all_2Dy,elliptic_2D,Ax,Ay,...
  eigen_all_2Dx_sf,eigen_all_2Dy_sf,elliptic_2D_sf,Ax_sf,Ay_sf,...
  eigen_all_SWx,eigen_all_SWy,elliptic_SW,Ax_SW,Ay_SW,...
  eigen_all_SWx_sf,eigen_all_SWy_sf,elliptic_SW_sf,Ax_SW_sf,Ay_SW_sf,...
  eigen_all_SWEx,eigen_all_SWEy,elliptic_SWE,Ax_SWE,Ay_SWE,...
  eigen_all_SWEx_sf,eigen_all_SWEy_sf,elliptic_SWE_sf,Ax_SWE_sf,Ay_SWE_sf...
  ]=ECT(flg,'gsd',gsd,'nodeState',nodeState,'sedTrans',sedTrans,'hiding',hiding,'ad',u_b,'secflow',[I,1,1],'twoD',twoDnr);
      
%% CFL

max_dt_qs=input.grd.dx/max(eigen_all_qs);
max_dt_fc=input.grd.dx/max(eigen_all);