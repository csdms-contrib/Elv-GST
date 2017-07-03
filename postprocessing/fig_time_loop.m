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
%$Id: fig_time_loop.m 106 2017-06-27 12:45:13Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/postprocessing/fig_time_loop.m $
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

function fig_time_loop(path_fold_main,fig_input)


%% 
%% READ
%% 

%paths
path_file_input=fullfile(path_fold_main,'input.mat');

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
            path_file_output_sng=fullfile(path_fold_temp_output,sprintf('%06d.mat',kT));
            output_par=load(path_file_output_sng);
            for ko=1:input.mdv.no
                aux_varname=input.mdv.output_var{1,ko}; %variable name to update in output.mat
                switch aux_varname
                    case 'time_loop'
                        output_m.(aux_varname)((kT-2)*floor(input.mdv.Flmap_dt/input.mdv.dt)+1:(kT-1)*floor(input.mdv.Flmap_dt/input.mdv.dt))=output_par.(aux_varname);
                    otherwise 
%                         nel=size(output_m.(aux_varname)); %size of the variable in the .mat file
%                         output_m.(aux_varname)(1:nel(1),1:nel(2),1:nel(3),kT)=output_par.(aux_varname)(1:nel(1),1:nel(2),1:nel(3));
                end
            end
        end
end

%% RENAME

%input
v2struct(input.mdv,{'fieldnames','time_results','nt','time'});

nT=numel(time_results);

%fig_input
v2struct(fig_input.tlo);

% aux_tl=reshape(output_m.time_loop,nt,nT);
% kT=find(isnan(aux_tl(1,2:end)));
% if isempty(kT) %get the last results
%     kT=nT;
% end
% time_loop_p=aux_tl(:,kT);
time_loop_p=output_m.time_loop;

%% FIGURE INITIALIZE

han.fig=figure('name',prnt.filename);
set(han.fig,'units','centimeters','paperposition',prnt.size)
% set(han.fig,'units','normalized','outerposition',[0,0,1,1]) %full monitor 1
% set(han.fig,'units','normalized','outerposition',[-1,0,1,1]) %full monitor 2
[mt,mb,mr,ml,sh,sv]=pre_subaxis(han.fig,marg.mt,marg.mb,marg.mr,marg.ml,marg.sh,marg.sv);

%subplots initialize
    %if regular
for kr=1:npr
    for kc=1:npc
        han.sfig(kr,kc)=subaxis(npr,npc,kc,kr,1,1,'mt',mt,'mb',mb,'mr',mr,'ml',ml,'sv',sv,'sh',sh);
    end
end
    %if irregular
% han.sfig(1,1)=subaxis(npr,npc,1,1,1,1,'mt',mt,'mb',mb,'mr',mr,'ml',ml,'sv',sv,'sh',sh);

    %add axis to sub1
% pos.sfig(1,1)=han.sfig1.Position; % position of first axes    
% han.sfig(1,2)=axes('Position',pos.sfig1,'XAxisLocation','top','YAxisLocation','right','Color','none');

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
han.sfig(kr,kc).Title.String=titlestr{kr,kc};
% han.sfig(kr,kc).XColor='r';
% han.sfig(kr,kc).YColor='k';

%% PLOTS

kr=1; kc=1;    
han.p(kr,kc,1)=plot(time_loop_p,'parent',han.sfig(kr,kc),'color','k','linewidth',prop.lw1,'linestyle',prop.ls1);

end %function



