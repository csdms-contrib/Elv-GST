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
%$Id: fig_level.m 79 2017-04-20 11:07:54Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/postprocessing/fig_level.m $
%
%function_name does this and that


%INPUT:
%   -
%
%OUTPUT:
%   -
%
%HISTORY:
%160223
%   -V. Created for the first time.

function fig_level(path_fold_main,fig_input)

%% 
%% READ
%% 

%paths
path_file_input=fullfile(path_fold_main,'input.mat');
% path_file_output=fullfile(path_fold_main,'output.mat');

%input (input)
input=NaN;
load(path_file_input); 

%output (output)
switch fig_input.mdv.wh
    case 1
        output_m=matfile(fullfile(path_fold_main,'output.mat')); %matfile io object creation
    case 2
        %load the empty results
        output_m=load(fullfile(path_fold_main,'output.mat'));

        %load the separate resutls files and copy to the variable with all the results
        path_fold_temp_output=fullfile(path_fold_main,'TMP_output');
        dir_temp_output=dir(path_fold_temp_output);
        nto=numel(dir_temp_output)-2;
        for kT=2:nto
            path_file_output_sng=fullfile(path_fold_temp_output,sprintf('%06d.mat',kT)); %path to the separate file with results
            output_par=load(path_file_output_sng); %load the partial results
            for ko=1:input.mdv.no %loop on varaibles to save
                aux_varname=input.mdv.output_var{1,ko}; %variable name to update in output.mat
                switch aux_varname 
%                     case 'time_loop' %history variable
%                         if kT==nto %only save the last one
%                             output_m.(aux_varname)=output_par.(aux_varname);
%                         end
                    case {'etab','h'}
                        nel=size(output_m.(aux_varname)); %size of the variable in the .mat file
                        output_m.(aux_varname)(1:nel(1),1:nel(2),1:nel(3),kT)=output_par.(aux_varname)(1:nel(1),1:nel(2),1:nel(3));
%                     otherwise %map variable

                end
            end
            fprintf('%4.1f files readed \n',kT/nto*100)
        end
end

%% RENAME

%input
v2struct(input.mdv,{'fieldnames','nx','xcen','time_results'});
% v2struct(input.sed,{'fieldnames','dk'});

nT=numel(time_results);

%due to the way some functions were prepared for D3D
% in.nx=nx+2; 

%fig_input
v2struct(fig_input.lev);

%time vector
switch time_input
    case 0 
        time_v=1:1:nT;
    case 1
        time_v=time;
    case 2
        time_v=1:time:nT;
end

%figure position
switch disppos
	case 1
		disppos_vec=[0,0,1,1];
	case 2
		disppos_vec=[-1,0,1,1];
	case 3
		disppos_vec=[-1,0,2,1];
end

%% FIGURE INITIALIZE

han.fig=figure('name',prnt.filename);
set(han.fig,'units','centimeters','paperposition',prnt.size)
set(han.fig,'units','normalized','outerposition',disppos_vec)
[mt,mb,mr,ml,sh,sv]=pre_subaxis(han.fig,marg.mt,marg.mb,marg.mr,marg.ml,marg.sh,marg.sv);

%subplots initialize
    %if regular
for kr=1:npr
    for kc=1:npc
        han.sfig(kr,kc)=subaxis(npr,npc,kc,kr,1,1,'mt',mt,'mb',mb,'mr',mr,'ml',ml,'sv',sv,'sh',sh);
    end
end

%properties
    %sub11
kr=1; kc=1;   
hold(han.sfig(kr,kc),'on')
% axis(han.sfig(kr,kc),'equal')
han.sfig(kr,kc).Box='on';
if lims.auto==0
    han.sfig(kr,kc).XLim=lims.x(kr,kc,:);
    han.sfig(kr,kc).YLim=lims.y(kr,kc,:);
end
han.sfig(kr,kc).XLabel.String=xlabels{kr,kc};
han.sfig(kr,kc).YLabel.String=ylabels{kr,kc};
% han.sfig(kr,kc).XTickLabel='';
% han.sfig(kr,kc).YTickLabel='';
% han.sfig(kr,kc).XTick=[];  
% han.sfig(kr,kc).YTick=[];  
% han.sfig(kr,kc).XScale='log';
% han.sfig(kr,kc).YScale='log';
han.sfig(kr,kc).Title.String=sprintf('%s=%f',titlestr{kr,kc},0);
% han.sfig(kr,kc).XColor='r';
% han.sfig(kr,kc).YColor='k';

%% FIRST UPDATE

kt=1;

%copy paste this section in the time loop
    %% X DATA

x=xcen*unitx;

    %% Y DATA

etab=output_m.etab(:,:,:,kt)*unity;
etaw=etab+output_m.h(:,:,:,kt)*unity;

%% PLOTS

kr=1; kc=1;    
han.p(kr,kc,1)=plot(x,etab,'parent',han.sfig(kr,kc),'color','k','linewidth',prop.lw1,'linestyle',prop.ls1);
han.p(kr,kc,2)=plot(x,etaw,'parent',han.sfig(kr,kc),'color','b','linewidth',prop.lw1,'linestyle',prop.ls1);

%% UPDATE FIGURE 

for kt=time_v

    %% Y DATA

etab=output_m.etab(:,:,:,kt)*unity;
etaw=etab+output_m.h(:,:,:,kt)*unity;

%% UPDATE FIGURE DATA

kr=1; kc=1;
han.p(kr,kc,1).YData=etab;
han.p(kr,kc,2).YData=etaw;

%% TITLE
han.sfig(kr,kc).Title.String=sprintf('%s=%3.2f',titlestr{kr,kc},time_results(kt)*unitt);

%% DISPLAY
switch dispfig
    case 1
        pause
    case 2
        pause(pausetime)
end
   

end

end %function



