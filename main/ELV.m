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
%$Id: ELV.m 108 2017-06-27 14:46:40Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/ELV.m $
%
%ELV is the main function of the model
%
%ELV(path_file_input,fid_log)
%
%INPUT:
%   -path_file_input = path to the file input.mat; [char]; 
%	-fid_log = log file identifier
%
%OUTPUT:
%   -
%
%HISTORY:
%160223
%   -V. Created for the first time.
%
%160415
%   -V. First Initial Condition and then Boundary Condition in case you
%   want equilibrium BC.
%
%160614
%   -V. Output pmm for adjusting CFL condition
%
%160429
%   -L. Introduce funtion conditions_construction
%
%160803
%	-L merged Vv3, Lv3
%
%161010
%   -L&V. Nourishment crap
%
%170126
%   -L. Repetitive nourishments
%
%170317
%   -V. Mean loop time
%
%170324
%   -V. Friction correction

function ELV(path_file_input,fid_log)
% version='5';
% fprintf(fid_log,'ELV version: %s\n',version);

%% 
%% INITIALIZATION
%% 

%% INPUT READING

fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of input reading');
input=NaN; %V is stupid and has just realised that 'input' is also a function in MatLab. GEFELICITEERD!
load(path_file_input); %(input)

%% INPUT CHECK

fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of input checking');
input=check_input(input,path_file_input,fid_log); 

%% 
%% PREPROCESSING
%% 

%% PREALLOCATE OTHER VARIABLES AND COUNTERS

nourcount = 1; %nourishment counter
ktl=1; %start time_loop counter
time_loop=zeros(floor(input.mdv.Flmap_dt/input.mdv.dt),1); %this needs to be before output_creation if you want to save the variable. Preallocate with zeros to avoid NaN in first saving loop.
tic_display=tic; %tic to control display in screen

%% INITIAL AND BOUNDARY CONDITION CONSTRUCTION

fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of initial and boundary condition construction');
[u,h,etab,Mak,La,msk,Ls,Cf,bc]=condition_construction(input,fid_log);

%% WRITE t0

fprintf(fid_log,'%s %s\n',datestr(datetime('now')),'Start of t0 conditions saving');
output_creation(input,fid_log)
kts=2; %saving time counter (in 1 there is the initial condition)

%% 
%% TIME LOOP
%% 

for kt=input.mdv.t0:input.mdv.nt   

tic_looptime=tic; %tic to track time spent in loop

%% FLOW UPDATE
[u,h]=flow_update(u,h,etab,Cf,bc,input,fid_log,kt);

%% FRICTION CORRECTION

[Cf_b]=friction_correction(u,h,Cf,Mak,La,input,fid_log,kt);

%% SEDIMENT TRANSPORT
[qbk,~]=sediment_transport(input.aux.flg,input.aux.cnt,h',(u(1,:).*h)',Cf_b,La',Mak',input.sed.dk,input.tra.param,input.aux.flg.hiding_parameter,1,fid_log,kt);qbk=qbk';

%% BED LEVEL UPDATE

etab_old=etab; %for Hirano
etab=bed_level_update(etab,qbk,bc,input,fid_log,kt);

%% ACTIVE LAYER THICKNESS UPDATE

La_old=La; %for Hirano
La=active_layer_thickness_update(h,Mak,La,input,fid_log,kt);

%% GRAIN SIZE DISTRIBUTION UPDATE

%save for the check
Mak_old=Mak; 
msk_old=msk;
Ls_old=Ls;

[Mak,msk,Ls,La,etab,ell_idx,celerities,pmm]=grain_size_distribution_update(Mak,msk,Ls,La_old,La,etab_old,etab,qbk,bc,u,h,Cf_b,input,fid_log,kt);

%% FRICTION UPDATE

Cf=friction(h,Mak,Cf,input,fid_log,kt);

%% CHECK SIMULATION

check_simulation(u,h,Mak,Mak_old,msk,msk_old,La,La_old,Ls,Ls_old,qbk,bc,ell_idx,celerities,pmm,input,fid_log,kt);

%% RESULTS WRITING
u = u(1,:);
if kts<=input.mdv.nT && kt*input.mdv.dt>=input.mdv.time_results(kts) %kt*dt>time_results(kts)
    print_tloop(time_loop,input,fid_log,kt,kts)
    write_results(input,fid_log,kts)
    kts=kts+1; %time save counter
    ktl=1; %reset time_loop counter
end

%% NOURISHMENT?
if kt*input.mdv.dt == input.nour.t(nourcount)
    [Mak,msk,Ls,La,etab]=add_nourishment(Mak,msk,Ls,La,etab,h,input,fid_log,kt);
    nourcount = nourcount + 1;
end

%% DISPLAY

time_loop(ktl)=toc(tic_looptime);
ktl=ktl+1; %update time_loop counter
if toc(tic_display)>input.mdv.disp_time
    display_tloop(time_loop,input,fid_log,kt,kts)
    tic_display=tic; %reset display tic
end

end %time loop

%% PUT OUTPUT FILES TOGETHER

join_results(input,fid_log);
