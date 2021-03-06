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
%$Id: fig_patch.m 106 2017-06-27 12:45:13Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/postprocessing/fig_patch.m $
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

function fig_patch(path_fold_main,fig_input)

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
v2struct(input.mdv,{'fieldnames','nx','nf','nsl','xedg','xcen','time_results'});
v2struct(input.sed,{'fieldnames','dk'});

nT=numel(time_results);

%due to the way some functions were prepared for D3D
in.nx=nx+2; 
in.nl=nsl+1; 

%fig_input
v2struct(fig_input.pat);

%time vector
switch time_input
    case 0 
        time_v=1:1:nT;
    case 1
        if isnan(time)
            time_v=nT;
        else
            time_v=time;
        end
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

%movie initialize
if fig_input.pat.dispfig==5
    paths_print=fullfile(path_fold_main,'figures',sprintf('%s',input.run));
    obj_vid=VideoWriter(paths_print);
    obj_vid.FrameRate=fig_input.pat.FrameRate;
    obj_vid.FrameRate=fig_input.pat.Quality;
    open(obj_vid);
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
han.sfig(kr,kc).Title.String=sprintf('%s=%f',titlestr{kr,kc},0);
% han.sfig(kr,kc).XColor='r';
% han.sfig(kr,kc).YColor='k';

%% FIRST UPDATE

kt=1;

%copy paste this section in the time loop
    %% X DATA

in.XCOR=xedg*unitx;

    %% Y DATA

% elevation of the substrate in.sub=[etab;etab-La;etab-La-Ls];
in.sub=[output_m.etab(:,:,:,kt);output_m.etab(:,:,:,kt)-output_m.La(:,:,:,kt);repmat(output_m.etab(:,:,:,kt),input.mdv.nsl,1)-repmat(output_m.La(:,:,:,kt),input.mdv.nsl,1)-reshape(cumsum(output_m.Ls(:,:,:,kt),3),input.mdv.nx,input.mdv.nsl,1,1)']*unity;

    %% C DATA

%put active layer and substrate in one matrix
M_all=NaN(nf,nx,nsl+1); %active layer and substrate in one matrix
M_all(1:nf-1,:,1    )=output_m.Mak(:,:,:,kt);
M_all(1:nf-1,:,2:end)=output_m.msk(:,:,:,kt);

L_all=NaN(nf,nx,nsl+1);
L_all(:,:,1      )=repmat(output_m.La(:,:,:,kt),nf,1,1,1);
L_all(:,:,2:nsl+1)=repmat(output_m.Ls(:,:,:,kt),nf,1,1,1);

F_all=M_all./L_all;
F_all(end,:,:)=ones(1,nx,nsl+1)-sum(F_all(1:nf-1,:,:),1);

% matrix to plot
switch cvar
    case 1 %volume fraction kf
        in.cvar=reshape(F_all(kf,:,:,:),nx,nsl+1)';
    case 2 %Dm        
        in.cvar=reshape(2.^(sum(F_all.*log2(repmat(dk,1,nx,nsl+1)/0.001),1)),nx,nsl+1)';
end 

    %% faces, vertices, color
[f,v,col]=rework4patch(in);

%% elliptic nodes
if input.mor.gsdupdate>1 || input.mor.ellcheck==1
    x_ell_idx=xcen(output_m.ell_idx(:,:,:,kt));
    y_ell_idx=aux_yelev(output_m.ell_idx(:,:,:,kt));
end

%% PLOTS

kr=1; kc=1;    
han.p(kr,kc,1)=patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'parent',han.sfig(kr,kc),'FaceColor','flat','EdgeColor',prop.edgecolor,'linewidth',prop.lw_patch);
if input.mor.gsdupdate>1 || input.mor.ellcheck==1
    han.s(kr,kc)=scatter(x_ell_idx,y_ell_idx,10,'m','filled');
end

%colormap
kr=1; kc=1;
view(han.sfig(kr,kc),[0,90]);
if lims.auto==0
    colormap(han.sfig(kr,kc),cmap);
    caxis(han.sfig(kr,kc),lims.c(kr,kc,1:2));
end

%colorbar
pos.sfig=han.sfig(cbar.sfig(1),cbar.sfig(2)).Position;
han.cbar=colorbar(han.sfig(cbar.sfig(1),cbar.sfig(2)),'location',cbar.location);
pos.cbar=han.cbar.Position;
han.cbar.Position=pos.cbar+cbar.displacement;
han.sfig(cbar.sfig(1),cbar.sfig(2)).Position=pos.sfig;
han.cbar.Label.String=cbar.label;

%% UPDATE FIGURE 

for kt=time_v
    %% X DATA

in.XCOR=xedg*unitx;

    %% Y DATA

% elevation of the substrate in.sub=[etab;etab-La;etab-La-Ls];
in.sub=[output_m.etab(:,:,:,kt);output_m.etab(:,:,:,kt)-output_m.La(:,:,:,kt);repmat(output_m.etab(:,:,:,kt),input.mdv.nsl,1)-repmat(output_m.La(:,:,:,kt),input.mdv.nsl,1)-reshape(cumsum(output_m.Ls(:,:,:,kt),3),input.mdv.nx,input.mdv.nsl,1,1)']*unity;

    %% C DATA

%put active layer and substrate in one matrix
M_all=NaN(nf,nx,nsl+1); %active layer and substrate in one matrix
M_all(1:nf-1,:,1    )=output_m.Mak(:,:,:,kt);
M_all(1:nf-1,:,2:end)=output_m.msk(:,:,:,kt);

L_all=NaN(nf,nx,nsl+1);
L_all(:,:,1      )=repmat(output_m.La(:,:,:,kt),nf,1,1,1);
L_all(:,:,2:nsl+1)=repmat(output_m.Ls(:,:,:,kt),nf,1,1,1);

F_all=M_all./L_all;
F_all(end,:,:)=ones(1,nx,nsl+1)-sum(F_all(1:nf-1,:,:),1);

% matrix to plot
switch cvar
    case 1 %volume fraction kf
        in.cvar=reshape(F_all(kf,:,:,:),nx,nsl+1)';
    case 2 %Dm        
        in.cvar=reshape(2.^(sum(F_all.*log2(repmat(dk,1,nx,nsl+1)/0.001),1)),nx,nsl+1)';
end 

    %% faces, vertices, color
[f,v,col]=rework4patch(in);

%% UPDATE FIGURE DATA

kr=1; kc=1;
han.p(kr,kc,1).Faces=f;
han.p(kr,kc,1).Vertices=v;
han.p(kr,kc,1).CData=col;

if input.mor.gsdupdate>1 || input.mor.ellcheck==1
    x_ell_idx=xcen(output_m.ell_idx(:,:,:,kt));
    y_ell_idx=aux_yelev(output_m.ell_idx(:,:,:,kt));
    
    han.s(kr,kc).XData=x_ell_idx;
    han.s(kr,kc).YData=y_ell_idx;
end

%% TITLE
han.sfig(kr,kc).Title.String=sprintf('%s=%3.2f',titlestr{kr,kc},time_results(kt)*unitt);


%% general
set(findall(han.fig,'-property','FontSize'),'FontSize',prop.fs)
% set(findall(han.fig,'-property','FontName'),'FontName',prop.fn) %!!! attention, there is a bug in Matlab and this is not enforced. It is necessary to change it in the '.eps' to 'ArialMT' (check in a .pdf)

%% DISPLAY 

switch dispfig
    case 1
        pause
    case 2
        pause(pausetime)
    case 3
        paths_print=fullfile(path_fold_main,'figures',sprintf('%s_%04d.png',prnt.filename,kt));
        print(han.fig,paths_print,'-dpng','-r300')
    case 4
        paths_print=fullfile(path_fold_main,'figures',sprintf('%s_%04d.eps',prnt.filename,kt));
        print(han.fig,paths_print,'-depsc2','-loose','-cmyk')
    case 5
        frame=getframe(han.fig);
        writeVideo(obj_vid,frame);
end
   
end

%movie close
if fig_input.pat.dispfig==5
    close(obj_vid);
end 

end %function



