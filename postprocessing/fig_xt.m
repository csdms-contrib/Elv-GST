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
%$Id: fig_xt.m 106 2017-06-27 12:45:13Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/postprocessing/fig_xt.m $
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
%
%160603
%   -V. Print .eps added

function fig_xt(path_fold_main,fig_input)

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
                    case {'etab','La','Mak','msk','Ls','ell_idx'}
                        nel=size(output_m.(aux_varname)); %size of the variable in the .mat file
                        output_m.(aux_varname)(1:nel(1),1:nel(2),1:nel(3),kT)=output_par.(aux_varname)(1:nel(1),1:nel(2),1:nel(3));
%                     otherwise %map variable

                end
            end
        end
end

%% RENAME

%input
v2struct(input.mdv,{'fieldnames','nx','nf','nsl','nT','xedg','xcen','time_results'});
v2struct(input.sed,{'fieldnames','dk'});

%fig_input
v2struct(fig_input.xtv);

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
    aux_yelev=lims.y(kr,kc,1)*ones(1,nx);
else
    min_elev=-0.15;
    aux_yelev=min_elev*ones(1,nx);
end
han.sfig(kr,kc).XLabel.String=xlabels{kr,kc};
han.sfig(kr,kc).YLabel.String=ylabels{kr,kc};
% han.sfig(kr,kc).XTickLabel='';
% han.sfig(kr,kc).YTickLabel='';
% han.sfig(kr,kc).XTick=[];  
% han.sfig(kr,kc).YTick=[];  
% han.sfig(kr,kc).XScale='log';
% han.sfig(kr,kc).YScale='log';
% han.sfig(kr,kc).Title.String=sprintf('%s=%f',titlestr{kr,kc},0);
% han.sfig(kr,kc).XColor='r';
% han.sfig(kr,kc).YColor='k';


    %% C DATA
switch varp
    case 1
        error('check, I have not tried!')
        %put active layer and substrate in one matrix
        M_all=NaN(nf,nx,nsl+1,nT); %active layer and substrate in one matrix
        M_all(1:nf-1,:,1    ,:)=output_m.Mak(:,:,:,:);
        M_all(1:nf-1,:,2:end,:)=output_m.msk(:,:,:,:);

        L_all=NaN(nf,nx,nsl+1,nT);
        L_all(:,:,1      ,:)=repmat(output_m.La(:,:,:,:),nf,1,1,1);
        L_all(:,:,2:nsl+1,:)=repmat(output_m.Ls(:,:,:,:),nf,1,1,1); %size(repmat(output_m.Ls(:,:,:,:),nf,1,1,1)) size(L_all(:,:,2:nsl+1,nT))

        F_all=M_all./L_all;
        F_all(end,:,:,:)=ones(1,nx,nsl+1,nT)-sum(F_all(1:nf-1,:,:,:),1);
        
        %get the one we want
        zp=squeeze(F_all(varp2,kxp,1,:));
    case 2
        error('check, I have not tried!')
        zp=squeeze(sum(output_m.qbk(:,kxp,:,:),1)).*unitqb;
        
    case 3
        error('check, I have not tried!')
        %qbk
    case 4
        zp=squeeze(output_m.pmm(varp2,:,:,:));
end


% matrix to plot
% switch cvar
%     case 1 %volume fraction kf
%         in.cvar=reshape(F_all(kf,:,:,:),nx,nsl+1)';
%     case 2 %Dm        
%         in.cvar=reshape(2.^(sum(F_all.*log2(repmat(dk,1,nx,nsl+1)/0.001),1)),nx,nsl+1)';
% end 

%% elliptic nodes
% if input.mor.gsdupdate>1 || input.mor.ellcheck==1
%     x_ell_idx=xcen(output_m.ell_idx(:,:,:,kt));
%     y_ell_idx=aux_yelev(output_m.ell_idx(:,:,:,kt));
% end

%% PLOTS

[xp,yp]=meshgrid(xcen.*unitx,time_results.*unitt);

kr=1; kc=1;    
han.p(kr,kc,1)=surf(xp,yp,zp',zp','parent',han.sfig(kr,kc),'edgecolor',prop.edgecolor);

% if input.mor.gsdupdate>1 || input.mor.ellcheck==1
%     han.s(kr,kc)=scatter(x_ell_idx,y_ell_idx,10,'m','filled');
% end

%colormap
% kr=1; kc=1;
% view(han.sfig(kr,kc),[0,90]);
% if lims.auto==0
%     colormap(han.sfig(kr,kc),cmap);
%     caxis(han.sfig(kr,kc),lims.c(kr,kc,1:2));
% end

%colorbar
pos.sfig=han.sfig(cbar.sfig(1),cbar.sfig(2)).Position;
han.cbar=colorbar(han.sfig(cbar.sfig(1),cbar.sfig(2)),'location',cbar.location);
pos.cbar=han.cbar.Position;
han.cbar.Position=pos.cbar+cbar.displacement;
han.sfig(cbar.sfig(1),cbar.sfig(2)).Position=pos.sfig;
han.cbar.Label.String=cbar.label;


%% TITLE
% han.sfig(kr,kc).Title.String=titlestr{kr,kc};

%% DISPLAY 

switch dispfig
    case 1
        paths_print=fullfile(path_fold_main,'figures',sprintf('%s_%d.png',prnt.filename,varp));
        print(han.fig,paths_print,'-dpng','-r300')
    case 2
        paths_print=fullfile(path_fold_main,'figures',sprintf('%s_%d.eps',prnt.filename,varp));
        print(han.fig,paths_print,'-depsc2','-loose','-cmyk')
end
   
end




