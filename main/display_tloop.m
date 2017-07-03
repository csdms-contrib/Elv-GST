%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 81 $
%$Date: 2017-04-20 15:51:02 +0200 (Thu, 20 Apr 2017) $
%$Author: V $
%$Id: display_tloop.m 81 2017-04-20 13:51:02Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/display_tloop.m $
%
%display_tloop is a function that displays in command at the end of every time step
%
%display_tloop(time_loop,input,fid_log,kt)
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
%160524
%   -V. Addition of the model run_id
%
%170317
%   -V. Mean loop time

function display_tloop(time_loop,input,fid_log,kt,kts)
%comment out fot improved performance if the version is clear from github
% version='2';
% if kt==1; fprintf(fid_log,'display_tloop version: %s\n',version); end 

%% RENAME
nt=input.mdv.nt;
run_id=input.run;
disp_t_nt=input.mdv.disp_t_nt;

%%

if kts>2
    kt_mean_idx=1:1:disp_t_nt;
elseif kt>disp_t_nt
    kt_mean_idx=kt-disp_t_nt:kt-1;
else
    kt_mean_idx=1:kt-1;
end

expected_s=(nt-kt)*mean(time_loop(kt_mean_idx)); %expected seconds to finish;

fl_d=floor(expected_s/3600/24); %floor of the remaining days
fl_h=floor((expected_s-fl_d*3600*24)/3600); %floor of the remaining hours
fl_m=floor((expected_s-fl_d*3600*24-fl_h*3600)/60); %floor of the remaining minutes
fl_s=expected_s-fl_d*3600*24-fl_h*3600-fl_m*60; %remaining minutes

fprintf('run %s, %4.2f%% computed, %01dd %02dh %02dm %2.0fs to finish \n',run_id,kt/nt*100,fl_d,fl_h,fl_m,fl_s);

