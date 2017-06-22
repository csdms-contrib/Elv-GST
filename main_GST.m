%main track front GST ELV results

%%
%% INPUT
%%

%simulations
in.run={'067','068'};
in.serie={'H','H'};

dt_a=365*24*3600; %simple model time step [s]  (if NaN use Flmap_dt)
t0_a=50*365*24*3600; %simple model time 0 [s]
x0_a=5000; %position of the front at time 0 [m]

% dt_a=1*24*3600; %simple model time step [s]  (if NaN use Flmap_dt)
% t0_a=15*365*24*3600; %simple model time 0 [s]
% x0_a=250; %position of the front at time 0 [m]

%%
%% CALC
%%

%folder
in.fold_pwd=pwd;
in.fold_runs=fullfile(in.fold_pwd,'sim');

%param
qb2_thres=1e-7; %minimum transpor of gravel to not be considered 0
nun=3; 

%initial values for iteration
ks=1;
known(ks,1).Qbk0=[0.7,0.3]*1e-5;
known(ks,1).h_u=1.202646;
known(ks,1).slopeb_u=9.376435e-04;
known(ks,1).Fak_u=0.1053120;

known(ks,1).h_d=4.017933;
known(ks,1).slopeb_d=2.514449e-05;
known(ks,1).Fak_d=1.0;

ks=2; %same
known(ks,1)=known(1,1);

%% PREALLOCATE

ns=numel(in.run); %number of simulations

%numerical
front_cord=cell(ns,1);
slope=cell(ns,1);
Fa1=cell(ns,1);
qb1_front=cell(ns,1);
thetak=cell(ns,1);
u_st=cell(ns,1);
dm_x=cell(ns,1);
time_results_all=cell(ns,1);
c_d=cell(ns,1);
t_d=cell(ns,1);

%analytical
c_a=cell(ns,1);
time_results_anl=cell(ns,1);
x_a=cell(ns,1);
slopeb_anl=cell(ns,1);
F_u_anl=cell(ns,1);
thetak_anl=cell(ns,1);
h_anl=cell(ns,1);

%% LOOP on simulations

for ks=1:ns

%% READ

%paths
path_fold_main=fullfile(in.fold_runs,in.serie{ks},in.run{ks});
path_file_input=fullfile(path_fold_main,'input.mat');

%input (input)
input=NaN;
load(path_file_input); 

%output (output)
output_m=matfile(fullfile(path_fold_main,'output.mat')); %matfile io object creation

%% RENAME

%input
v2struct(input.mdv,{'fieldnames','nT','xcen','Flmap_dt','time_results','nx','nef'});
v2struct(input.grd,{'fieldnames','dx'});
v2struct(input.sed,{'fieldnames','dk'});
v2struct(input.mor,{'fieldnames','interfacetype','fIk_alpha'});

%output
qbk=output_m.qbk;
etab=output_m.etab;
Fak=output_m.Mak./repmat(output_m.La,nef,1,1,1);
h=output_m.h;
u=output_m.u;
Cf=output_m.Cf;

time_results_all{ks,1}=time_results;

%% ANALYSIS OF RUN

front_cord{ks,1}=NaN(nT-1,1);
for kT=2:nT
    front_cord_t=find(qbk(2,:,1,kT),1,'last'); %search for the last point where gravel is in transport
    if isempty(front_cord_t) || isnan(front_cord_t) %
        front_cord{ks,1}(kT,1)=1;
    elseif front_cord_t>=nx-1
        front_cord{ks,1}(kT,1)=find(qbk(2,:,1,kT)>qb2_thres,1,'last'); %search for the last point where gravel is over threshold
    else
        front_cord{ks,1}(kT,1)=front_cord_t;
    end

end %kT

%% celerity
[x_n,t_n,~]=unique(front_cord{ks,1});
c_n=diff(x_n(1:end-1))./diff(t_n(1:end-1)); %celerity  [space_nodes/time_nodes]
c_d{ks,1}=c_n.*dx./Flmap_dt; %dimensional celerity
t_d{ks,1}=t_n(1:end-2).*Flmap_dt;

%% derived data

% preallocate
slope{ks,1}=NaN(nT-1,2);
Fa1{ks,1}=NaN(nT-1,2);
qb1_front{ks,1}=NaN(nT-1,1);
thetak{ks,1}=NaN(nT-1,2,2);
u_st{ks,1}=NaN(nT-1,nx);
dm_x{ks,1}=NaN(nT-1,nx);

% loop on time
for kT=2:nT
    if front_cord{ks,1}(kT,1)-nun<=nun
        slope{ks,1}(kT,1)=NaN;
        Fa1{ks,1}(kT,1)=NaN;
        qb1_front{ks,1}(kT,1)=NaN;
        thetak{ks,1}(kT,:,:)=NaN;
        u_st{ks,1}(kT,1)=NaN;
    else
        %slope
        aux_slope_u=polyfit(xcen(1:front_cord{ks,1}(kT,1)-nun),etab(1,1:front_cord{ks,1}(kT,1)-nun,1,kT),1);
        aux_slope_d=polyfit(xcen(front_cord{ks,1}(kT,1)+nun:end),etab(1,(front_cord{ks,1}(kT,1)+nun:end),1,kT),1);
        slope{ks,1}(kT,:)=-[aux_slope_u(1),aux_slope_d(1)];
        %fractions
        Fak_x=Fak(1,:,1,kT);
        Fa_x=[Fak_x;1-sum(Fak_x,1)];
        Fak_x_u=Fak(1,1:front_cord{ks,1}(kT,1)-nun,1,kT);
        Fak_x_d=Fak(1,front_cord{ks,1}(kT,1)+nun:end,1,kT);
        Fa1_u=mean(Fak_x_u);
        Fa1_d=mean(Fak_x_d);
        Fa1{ks,1}(kT,:)=[Fa1_u,Fa1_d];
        %qb1
        qb1_front{ks,1}(kT,1)=qbk(1,front_cord{ks,1}(kT,1),1,kT);
        %shields
        thetak_u=mean(Cf(1,1:front_cord{ks,1}(kT,1)-nun,1,kT).*u(1,1:front_cord{ks,1}(kT,1)-nun,1,kT).^2./1.65./9.81)./dk; %att! v<2016b will not allow this operation
        thetak_d=mean(Cf(1,front_cord{ks,1}(kT,1)+nun:end,1,kT).*u(1,front_cord{ks,1}(kT,1)+nun:end,1,kT).^2./1.65./9.81)./dk; %att! v<2016b will not allow this operation
        thetak{ks,1}(kT,1,:)=thetak_u;
        thetak{ks,1}(kT,2,:)=thetak_d;
        %u star
        u_st{ks,1}(kT,:)=sqrt(Cf(1,:,1,kT)).*u(1,:,1,kT);
        %Dm
        dm_x{ks,1}(kT,:)=2.^(sum(Fa_x.*log2(repmat(dk,1,nx)),1)); 
    end %front
end %kT

%% SIMPLE MODEL

%Hoey94 in case of Hirano
switch interfacetype
    case 1
        fIk_alpha=0;
end

%rename input tra
rename_tra

if isnan(dt_a)
    dt_a=Flmap_dt;
end
time_results_anl{ks,1}=(t0_a:dt_a:time_results_all{ks,1}(end))';
nT_a=numel(time_results_anl{ks,1});
c_a{ks,1}=NaN(nT_a,1);
x_a{ks,1}=NaN(nT_a,1);
    x_a{ks,1}(1,1)=x0_a;
slopeb_anl{ks,1}=NaN(nT_a,2);
F_u_anl{ks,1}=NaN(nT_a,1);    
thetak_anl{ks,1}=NaN(nT_a,2,2);    

%boundary condition for each saved time
Qbk0=NaN(nT_a,input.mdv.nf);
for kf=1:input.mdv.nf
    Qbk0(:,kf)=interp1(input.bcm.timeQbk0,input.bcm.Qbk0(:,kf),time_results_anl{ks,1},'linear')'; 
end
bc.qbk0 = Qbk0./input.grd.B(1);

%loop on time
for kp=1:nT_a

time_l=time_results_anl{ks,1}(kp,1);
timeQbk0_a_coor=kp;

%upstream slope
input_u=input; %copy input of the upstream end
input_u.bcm.Qbk0=bc.qbk0(timeQbk0_a_coor,:); %full load
%sand reach slope
input_d=input; %copy input of the upstream end
Qb0=sum(bc.qbk0(timeQbk0_a_coor,:),2);
pg0=bc.qbk0(timeQbk0_a_coor,2)./Qb0;

%initial values for iteration
kb=1;
    %upstream
input_u.ini.Fak=known(ks,kb).Fak_u; %effective fractions at the active layer [-]; [(nf-1)x1 double] | [(nf-1)xnx double]; e.g. [0.2;0.3]
input_u.ini.h=known(ks,kb).h_u; %flow depth [m]; [1x1 double] | [1xnx double]; e.g. [1.5]
input_u.ini.slopeb=known(ks,kb).slopeb_u; %bed slope [-]; [1x1 double] | [1xnx double]; e.g. [1e-3]
    %downstream
input_d.ini.Fak=known(ks,kb).Fak_d; %effective fractions at the active layer [-]; [(nf-1)x1 double] | [(nf-1)xnx double]; e.g. [0.2;0.3]
input_d.ini.h=known(ks,kb).h_d; %flow depth [m]; [1x1 double] | [1xnx double]; e.g. [1.5]
input_d.ini.slopeb=known(ks,kb).slopeb_d; %bed slope [-]; [1x1 double] | [1xnx double]; e.g. [1e-3]

%upstream
[slopeb_u,Fa_u,sum_Fak_obj,max_rel_error]=get_equivals_V(input_u.bch.Q0(1),input_u,sum(input_u.bcm.Qbk0),input_u.bcm.Qbk0/sum(input_u.bcm.Qbk0));
if max_rel_error>0.001 || sum_Fak_obj>1e-8
   warning('problem')
end
Fak_u=Fa_u(1);
h_u=(input_u.mdv.Cf.*(input_u.bch.Q0(1)/input_u.grd.B).^2/(9.81*slopeb_u)).^(1/3);
u_u=(input_u.bch.Q0(1)/input_u.grd.B)/h_u;
F_u=1-Fak_u; %surface gravel fraction at the upstream reach
thetak_u=(repmat(h_u*slopeb_u/1.65,2,1)./input.sed.dk)';

%downstream
Qb1_gst=Qb0*(1-pg0/F_u); %sand at the GST [m^3/s]
input_d.bcm.Qbk0=[Qb1_gst,0];
[slopeb_d,Fa_d,sum_Fak_obj,max_rel_error]=get_equivals_V(input_d.bch.Q0(1),input_d,sum(input_d.bcm.Qbk0),input_d.bcm.Qbk0/sum(input_d.bcm.Qbk0));
if max_rel_error>0.001 || sum_Fak_obj>1e-8
   warning('problem')
end
Fak_d=Fa_d(1);
h_d=(input_d.mdv.Cf.*(input_d.bch.Q0(1)/input_d.grd.B).^2/(9.81*slopeb_d)).^(1/3);
u_d=(input_d.bch.Q0(1)/input_d.grd.B)/h_d;
thetak_d=(repmat(h_d*slopeb_d/1.65,2,1)./input.sed.dk)';
    

F_w=fIk_alpha*pg0+(1-fIk_alpha)*F_u;

%base level 
time_etaw0_a=find(input.bch.timeetaw0<time_l,1,'last');

xi=(input.bch.etaw0(time_etaw0_a+1)-input.bch.etaw0(time_etaw0_a))/(input.bch.timeetaw0(time_etaw0_a+1)-input.bch.timeetaw0(time_etaw0_a));

if kp==1
    dh_u_dt=0;
    dslopeb_d_dt=0;
    dslopeb_u_dt=0;
else
    dh_u_dt=(h_u-h_anl{ks,1}(kp-1,1))/dt_a;
    dslopeb_d_dt=(slopeb_d-slopeb_anl{ks,1}(kp-1,2))/dt_a;
    dslopeb_u_dt=(slopeb_u-slopeb_anl{ks,1}(kp-1,1))/dt_a;
end

%celerity
c_a{ks,1}(kp)=1/(slopeb_u-slopeb_d)*(bc.qbk0(timeQbk0_a_coor,2)/(1-input.mor.porosity)/F_w/input.grd.B/x_a{ks,1}(kp)-xi+dh_u_dt-(input.grd.L-x_a{ks,1}(kp))*dslopeb_d_dt-0.5*x_a{ks,1}(kp)*dslopeb_u_dt); %celerity*distance 

%update position
x_a{ks,1}(kp+1)=x_a{ks,1}(kp)+c_a{ks,1}(kp)*dt_a;

aux.str=sprintf('%4.1f%% done of simulation %d of %d',kp/nT_a*100,ks,ns);
disp(aux.str);

%save
h_anl{ks,1}(kp,:)=[h_u,h_d];
slopeb_anl{ks,1}(kp,:)=[slopeb_u,slopeb_d];
F_u_anl{ks,1}(kp)=F_u;
thetak_anl{ks,1}(kp,1,:)=reshape(thetak_u,1,1,2);
thetak_anl{ks,1}(kp,2,:)=reshape(thetak_d,1,1,2);

end %kp
end %ks


%%
%% PLOT
%%
       
%% GST prop

fact.time=1/3600/24/365;

prop.color=[... %<  matlab 2014b default
 0.0000    0.0000    1.0000;... %blue
 0.0000    0.5000    0.0000;... %green
 1.0000    0.0000    0.0000;... %red
 0.0000    0.7500    0.7500;... %cyan
 0.7500    0.0000    0.7500;... %purple
 0.7500    0.7500    0.0000;... %ocre
 0.2500    0.2500    0.2500];   %grey

nr=1; nc=4;
figure('units','centimeters','paperposition',[0,0,27,15])

subplot(nr,nc,1)
hold on
for ks=1:ns
    han.p(ks)=plot(time_results_all{ks,1}.*fact.time,front_cord{ks,1}*dx,'color',prop.color(ks,:));
    han.p2(ks)=plot(time_results_anl{ks,1}.*fact.time,x_a{ks,1}(2:end),'color',prop.color(ks,:),'linestyle','--');
end
xlabel('time [years]')
ylabel('gst x-position [m]')

h=subplot(nr,nc,2);
hold on
for ks=1:ns
    plot(t_d{ks,1}.*fact.time,c_d{ks,1},'color',prop.color(ks,:));
    han.p2(ks)=plot(time_results_anl{ks,1}.*fact.time,c_a{ks,1},'color',prop.color(ks,:),'linestyle','--');
end
xlabel('time [years]')
ylabel('gst celerity [m/s]')

subplot(nr,nc,3)
hold on
for ks=1:ns
    plot(time_results_all{ks,1}.*fact.time,slope{ks,1}(:,1),'color',prop.color(ks,:))
    plot(time_results_all{ks,1}.*fact.time,slope{ks,1}(:,2),'color',prop.color(ks,:))
end
xlabel('time [years]')
ylabel('slope [-]')

subplot(nr,nc,4)
hold on
for ks=1:ns
    plot(time_results_all{ks,1}.*fact.time,1-Fa1{ks,1}(:,1),'color',prop.color(ks,:))
    plot(time_results_all{ks,1}.*fact.time,1-Fa1{ks,1}(:,2),'color',prop.color(ks,:))
end
xlabel('time [years]')
ylabel('Fg [-]')

