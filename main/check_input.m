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
%$Id: check_input.m 108 2017-06-27 14:46:40Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/check_input.m $
%
%check_input is a function that checks that the input is enough and makes sense
%
%input_out=check_input(input,path_file_input,fid_log)
%
%INPUT:
%   -input = variable containing the input [struct] e.g. input
%
%OUTPUT:
%   -input = variable containing the input [struct] e.g. input
%
%HISTORY:
%160223
%   -V. Created for the first time.
%
%160429
%   -L. changed hydrodynamic boundary conditions;
%       - new case numbers
%       - repeated hydrograph
%       - adjusted warning/error messages related to bc conditions
%       - boundary conditions from file
%160429
%   -L. adjust for morphodynamic boundary conditions;
%
%160429
%   -L adjust for more possibilities for initial conditions
% line 92 warning instead of error
%
%160705
%   -L entered path_folder_main check
%
%160818
%   -L put input.mor.interfacetype into nf>1 check
%   -L: line 468 typo?
%
%161010
%   -L&V. Nourishment crap
%   -L. case 13 1604011
%
%170123
%   -L. major updates based on updated manual


function input_out=check_input(input,path_file_input,fid_log)

%% 
%% RUN
%% 

if isfield(input,'run')==0 %if it does not exist 
    error('There is no run tag. Provide input.run')
else
    if ischar(input.run)==0
        error('input.run must be a char')
    end
end

%% 
%% MDV
%% 

input.mdv.nt=round(input.mdv.Tstop/input.mdv.dt); %number time steps [-] 
input.mdv.nx=round(input.grd.L/input.grd.dx); %number of grid points in space [-] 
input.mdv.nf=numel(input.sed.dk); %number of sediment size fraction [-]
input.mdv.nef=input.mdv.nf-1; %number of effective sediment size fraction [-]
input.mdv.no=numel(input.mdv.output_var); %number of output variables [-]

if input.mdv.nf~=1
    input.mdv.nsl=round(input.mor.total_ThUnLyr/input.mor.ThUnLyr); %number of vertical layers for substrate discretisation [-] 
else
    input.mdv.nef=1;
    input.mdv.nsl=1;
    input.mor.ThUnLyr=NaN;
    input.mor.La=NaN; %active layer thickness [m]; [1x1 double]; e.g. [0.1]
    input.ini.Fak=NaN; %effective fractions at the active layer [-]; [(nf-1)x1 double] | [(nf-1)xnx double]; e.g. [0.2;0.3]
    input.ini.fsk=NaN; %effective fractions at the substrate [-]; [(nf-1)x1 double] | [(nf-1)xnx double]; e.g. [0.2;0.3]
end

input.mdv.time=0:input.mdv.dt:input.mdv.nt*input.mdv.dt; %time vector [s]
input.mdv.time_results=0:input.mdv.Flmap_dt:input.mdv.nt*input.mdv.dt; %saving time vector [s]
input.mdv.xedg=0:input.grd.dx:input.mdv.nx*input.grd.dx; %cell edges vector x coordinate vector [m]
input.mdv.xcen=input.grd.dx/2:input.grd.dx:input.mdv.nx*input.grd.dx-input.grd.dx/2; %cell centers vector x coordinate vector [m]

input.mdv.nT=numel(input.mdv.time_results); %number of results [-]

% input.mdv.Tstop_round=nt*input.mdv.dt; %real computational time [s]
% input.mdv.L_round=nx*input.mdv.dx;%real domain length [m]

%scheme
if isfield(input.mdv,'UpwFac')==0
    input.mdv.UpwFac=1;
end
if input.mdv.UpwFac~=1
    warningprint(fid_log,'You are using a non-upwind scheme for Exner and Hirano. For a positive charactersitc celerity this is unstable, you better know what you are doing...')
end

%email
if isfield(input.mdv,'email')==0 %if it does not exist
    input.mdv.email.send=0; %0=NO; 1=YES
end
if input.mdv.email.send==1
    if isnumeric(input.mdv.email.password_coded)==0
        error('encode your password in input.mdv.email.password_coded')
    end
end

%save method
if isfield(input.mdv,'savemethod')==0 %if it does not exist
    input.mdv.savemethod=1; %1=directly in one file; 2=first in individual files
end

%starting time step
if isfield(input.mdv,'t0')==0 %if it does not exist
    input.mdv.t0=1; 
end

%display time 
if isfield(input.mdv,'disp_time')==0 %if it does not exist
    input.mdv.disp_time=1;
end
if isunix==1
    warningprint(fid_log,'It seems you are using a Unix machine, are you using the cluster? If the answer is "yes I do" the you should set the parameter input.mdv.disp_time equal to input.mdv.Tstop+1 so that there is no screen display')
end

%check
if isfield(input.mdv,'chk')==0 
    input.mdv.chk.mass=1; %mass check [-]; [1x1 double]; e.g. [1]
    input.mdv.chk.dM_lim=1e-5; %mass error limit [m^2, -!!]; [1x1 double]; e.g. [1e-8]
    input.mdv.chk.flow=1; %Froude and CFL check [-]; [1x1 double]; e.g. [1]
    input.mdv.chk.Fr_lim=0.8; %Fr limit [-]; [1x1 double]; e.g. [0.8]
    input.mdv.chk.cfl_lim=0.95; %CFL limit [-]; [1x1 double]; e.g. [0.95]
    input.mdv.chk.F_lim=0.01; %maximum error in volume fractions [-]; [1x1 double]; e.g. [0.01]
    input.mdv.chk.nan=1; %check for NaN in variables 0=NO; 1=YES;
    input.mdv.chk.ell=0; %display check for ellipticity 0=NO; 1=YES;
end

if isfield(input.mdv.chk,'mass')==0 
    input.mdv.chk.mass=1; %mass check [-]; [1x1 double]; e.g. [1]
    input.mdv.chk.dM_lim=1e-5; %mass error limit [m^2, -!!]; [1x1 double]; e.g. [1e-8]
else
    if isfield(input.mdv.chk,'dM_lim')==0 
        warningprint(fid_log,'You want to check if you lose mass but you are not specifying a limit (input.mdv.chk.dM_lim), I am using the default value of 1e-5')
        input.mdv.chk.dM_lim=1e-5; %mass error limit [m^2, -!!]; [1x1 double]; e.g. [1e-8]
    end
end

if isfield(input.mdv.chk,'flow')==0 
    input.mdv.chk.flow=1; %Froude and CFL check [-]; [1x1 double]; e.g. [1]
    input.mdv.chk.Fr_lim=0.8; %Fr limit [-]; [1x1 double]; e.g. [0.8]
    input.mdv.chk.cfl_lim=0.95; %CFL limit [-]; [1x1 double]; e.g. [0.95]
else
    if isfield(input.mdv.chk,'Fr_lim')==0 
        warningprint(fid_log,'You want to check the Froude number but you are not specifying a limit (input.mdv.chk.Fr_lim), I am using the default value 0.80')
        input.mdv.chk.Fr_lim=0.8; %Fr limit [-]; [1x1 double]; e.g. [0.8]
    end
    if isfield(input.mdv.chk,'cfl_lim')==0 
        warningprint(fid_log,'You want to check the CFL number but you are not specifying a limit (input.mdv.chk.cfl_lim), I am using the default value of 0.95')
        input.mdv.chk.cfl_lim=0.95; %CFL limit [-]; [1x1 double]; e.g. [0.95]
    end
end

if isfield(input.mdv.chk,'F_lim')==0 
    input.mdv.chk.F_lim=0.01; %maximum error in volume fractions [-]; [1x1 double]; e.g. [0.01]
end    
    
if isfield(input.mdv.chk,'nan')==0 
    input.mdv.chk.nan=1; %check for NaN in variables 0=NO; 1=YES;
end 

if isfield(input.mdv.chk,'ell')==0 
    input.mdv.chk.ell=0; %display check for ellipticity 0=NO; 1=YES;
end 

%% CONSTANTS

if isfield(input.mdv,'g')==0 
    input.mdv.g=9.81;
end
if isfield(input.mdv,'nu')==0 
    input.mdv.nu=1e-6;
end
if isfield(input.mdv,'rhow')==0 
    input.mdv.rhow=1000;
end
% if isfield(input.frc,'nk')==0 
%     input.frc.nk=2;
% end
if isfield(input.sed,'rhos')==0 
    input.sed.rhos=2650;
end

if isfield(input.mdv,'dd')==0 
    input.mdv.dd=1e-8; %diferential
end

%% PATHS

if isfield(input.mdv,'path_folder_main')==0 %if it does not exist the path to the main file
    [path_folder_main,~,~]=fileparts(path_file_input); %get path to main folder
    input.mdv.path_folder_main=path_folder_main;
end

if isfield(input.mdv,'path_folder_results')==0 %if it does not exist the path to the results folder
    input.mdv.path_folder_results=input.mdv.path_folder_main; %path to result folder 
end

if isfield(input.mdv,'path_file_output')==0 %if it does not exist the path to the output file
    input.mdv.path_file_output=fullfile(input.mdv.path_folder_main,'output.mat'); %path to the output file
end

%attention! this two folders are hardcoded in folders_creation because that function is called before check_input. Do not change (or do it with care :D )
%this should be changed by changint the function order
input.mdv.path_folder_TMP_output=fullfile(input.mdv.path_folder_main,'TMP_output'); %path to the temporal results output
input.mdv.path_folder_figures=fullfile(input.mdv.path_folder_main,'figures'); %path to the figures files

%% FRICTION

switch input.mdv.frictiontype
    case 1 %constant
        if isfield(input.mdv,'Cf')==0
            error('You need to specify input.mdv.Cf')
        end       
    otherwise
        error('input.mdv.frictiontype can be: 1 constant')
end

%wall correction
if isfield(input,'frc')==0
    input.frc.wall_corr=0; %flume wall correction: 0=NO; 1=YES
    input.frc.ripple_corr=0; %ripple correction: 0=NO; 1=YES
end

switch input.frc.wall_corr
    case 0
        input.frc.Cfb=input.mdv.Cf;
    case 1 
        warningprint(fid_log,'input.frc.wall_corr=1 has not been tested. Attention to the initial condition. Cfb needs to be passed to create it.');
        if isfield(input.frc,'H')==0 || isfield(input.frc,'L')==0
            error('You need to specify a ripple height (input.frc.H) and length (input.frc.L)')
        end
    case 2
        warningprint(fid_log,'input.frc.wall_corr=2 has not been tested. Attention to the initial condition. Cfb needs to be passed to create it.');
    case 3
        if isfield(input.frc,'Cfb')==0
           error('If you want to impose a value for the bed friction, you better input it...  input.frc.Cfb')
        end
end


switch input.frc.ripple_corr
    case 1 
        if isfield(input.frc,'H')==0 || isfield(input.frc,'L')==0
            error('You need to specify a ripple height (input.frc.H) and length (input.frc.L)')
        end
end
if input.mdv.flowtype == 6 %main-flood1-flood2
    if numel(input.mdv.Cf)==1
        if input.grd.crt==1
            warningprint(fid_log,'Only one friction coefficient is supplied, but the channel is rectangular.');
            warningprint(fid_log,'No action is taken.');
        else
            warningprint(fid_log,'Only one friction coefficient is supplied, it is used everywhere');
            input.mdv.Cf = [input.mdv.Cf; input.mdv.Cf; input.mdv.Cf];
        end
    elseif numel(input.mdv.Cf)==2
        warningprint(fid_log,'Only two friction coefficients are supplied, it is assumed you want both floodplains with same coefficient');
        input.mdv.Cf = [input.mdv.Cf input.mdv.Cf; input.mdv.Cf(2)];
    elseif numel(input.mdv.Cf)==3         
    else
        error('Wrong dimensions for input.mdv.Cf');
    end 
end

%%
%% WIDTH
%%

% Adjust input for the correct flowtype
switch input.mdv.flowtype
    case 6
        if isfield(input.grd,'crt')==0 %constant cross section?
            input.grd.crt=0; %not a constant cross section
        end
        
        if numel(input.grd.B)==1
            warningprint(fid_log,'You are selecting a constant width, it is beter switch to flowtype 1');
            input.mdv.flowtype = 1;
        end
        
        % Check for a main channel or also flood plaines
        if size(input.grd.B,1)==1
            if input.grd.crt==1 %constant rectangular cross section
                warningprint(fid_log,'We will use the variable width solver for rectangular channels');
                warningprint(fid_log,'The hydraulic radius is approximated as the flow depth');
                B = [input.grd.B,input.grd.B(end)]; 
                input.grd.dBdx = (B(2:end)-B(1:end-1))/(input.grd.dx); 
         
            else
                warningprint(fid_log,'We will use only the main channel, flood plane width is set to zero');
                input.grd.B = [input.grd.B; zeros(size(input.grd.B)); zeros(size(input.grd.B))];
                input.grd.Bparam = repmat([999666999; 0; 999666999; 0],1,size(input.grd.B,2));
            end
        elseif size(input.grd.B,1)==2
            warningprint(fid_log,'We will assume that there is only one flood plain');
            input.grd.B = [input.grd.B(1,:); input.grd.B(2,:); zeros(size(input.grd.B(1,:)))];
        elseif size(input.grd.B,1)==3
            % correct input
        else
            error('wrong input in width');
        end
      
        if input.grd.crt==1
        else
            if size(input.grd.Bparam,1)==2
                    input.grd.Bparam = [input.grd.Bparam; 999666999*ones(size(input.grd.Bparam(1,:))),zeros(size(input.grd.Bparam(1,:)))];   
            elseif size(input.grd.Bparam,1)==4
            else
                error('wrong input');
            end
        end
        
        % Interpolation part
        if size(input.grd.B,2)==1
            warningprint(fid_log,'Only one cross-section is specified, it is used everywhere');
            input.grd.B = repmat(input.grd.B,1,input.mdv.nx);
        elseif size(input.grd.B,2)<input.mdv.nx
            warningprint(fid_log,'Not enough cross-sections are specified, looking for an x-vector to apply interpolation');
            input.grd.B=interp1(input.grd.Bx,input.grd.B,repmat((1:input.mdv.nx)',3,1)); 
        elseif size(input.grd.B,2)>input.mdv.nx
            error('wrong specification of the width');
        end

        if input.grd.crt==1
        else
            if size(input.grd.Bparam,2)==1
                warningprint(fid_log,'Only one set of cross-section parameters is specified, it is used everywhere');
                input.grd.Bparam= repmat(input.grd.Bparam,1,input.mdv.nx);
            elseif size(input.grd.B,2)<input.mdv.nx
                warningprint(fid_log,'Not enough cross-sections parameters are specified, looking for an x-vector to apply interpolation');
                if isfield(input.grd,'Bxparam')
                elseif isfield(input.grd,'Bx')
                   input.grd.Bxparam = input.grd.Bx;
                else
                    error('Not enough input is specified');
                end           
                input.grd.Bparam = [interp1(input.grd.Bxparam,input.grd.Bparam(1,:),1:input.mdv.nx); interp1(input.grd.Bxparam,input.grd.Bparam(2,:),1:input.mdv.nx); interp1(input.grd.Bxparam,input.grd.Bparam(3,:),1:input.mdv.nx); ];       
            elseif size(input.grd.B,2)>input.mdv.nx
                error('wrong specification of the width');
            end
        end   
    otherwise %flow type not 6
        if size(input.grd.B,2)==1
            input.grd.B=repmat(input.grd.B,1,input.mdv.nx);
        else
            if size(input.grd.B,2)~=input.mdv.nx
                error('ups... seems your input is wrong');
            end
        end
end

%%
%% FOR FLOWTYPE 6, with variable width; check input
%% MDV and GRD
%%

switch input.mdv.flowtype
    case {1,5}
		if input.ini.initype~=1
			input.ini.u=NaN;
			input.ini.h=NaN;
		end
    case {2,3,4}
        if isfield(input.ini,'u')==0 || isfield(input.ini,'h')==0
            warningprint(fid_log,'For quasi-steady or unsteady flow type you must specify inital u and h')
        end
    case 6  
        if input.ini.initype~=1
			input.ini.u=NaN;
			input.ini.h=NaN;
        end       
end
    
%% DISPLAY TIME

if isfield(input.mdv,'disp_t_nt')==0
    input.mdv.disp_t_nt=floor(input.mdv.Flmap_dt/input.mdv.dt);
end
if mod(input.mdv.disp_t_nt,1)~=0
    warningprint(fid_log,'The value you input in input.mdv.disp_t_nt is not an integer. I have rounded it.')
    input.mdv.disp_t_nt=round(input.mdv.disp_t_nt);
end
if input.mdv.disp_t_nt>floor(input.mdv.Flmap_dt/input.mdv.dt)
    warningprint(fid_log,'The value you ask to average the time needed in each loop is larger than the saving time, I do not like you to waste computational time. I set it to the maximum. Check input.mdv.disp_t_nt')
    input.mdv.disp_t_nt=floor(input.mdv.Flmap_dt/input.mdv.dt);
end


%% 
%% INI
%%

%check that the kind of input is specified
if isfield(input.ini,'initype')==0
    error('You need to specify the kind of initial condition in input.ini.initype')
end

%check that the input matches the kind of input
switch input.ini.initype
    case 1 %normal flow (for a given qbk0)
        if isfield(input.bch,'timeQ0')==0 
            if input.bch.uptype == 12 && input.bch.dotype == 12
            else
                error('You need to provide the time at which the specific water discharge is specified (input.bch.timeQ0) if you want such an upstream hydrodynamic boundary condition')
            end
        elseif isfield(input.bch,'Q0')==0
            if input.bch.uptype == 12 && input.bch.dotype == 12
            else
                error('You need to provide the specific water discharge at the specified times (input.bch.Q0) if you want such an upstream hydrodynamic boundary condition')
            end
        end
        if isfield(input.bcm,'timeQbk0')==0
            if input.bcm.type == 12
            else
                error('You need to provide the time at which the volume of sediment transported excluding pores per unit time, and per size fraction is specified (input.bcm.timeQbk0) if you want such a morphodynamic boundary condition')
            end
        elseif isfield(input.bcm,'Qbk0')==0
            if input.bcm.type == 12
            else
                error('You need to provide the volume of sediment transported excluding pores per unit time, and per size fraction at the specified times (input.bcm.Qbk0) if you want such a morphodynamic boundary condition')
            end
        end
%         if isfield(input.ini,'Fak')==0 || isfield(input.ini,'u')==0 || isfield(input.ini,'h')==0 || isfield(input.ini,'slopeb')==0
%             error('You need to specify an initial value of Fak, u, h, and slopeb as an initial guess to find the normal flow intial condition')
%         end


    case 2 %free
        %check dimension if vectors
        if     numel(input.ini.u)~=1 && size(input.ini.u,2)~=input.mdv.nx
            error('The dimensions of input.ini.u are incorrect')
        elseif numel(input.ini.h)~=1 && size(input.ini.h,2)~=input.mdv.nx
            error('The dimensions of input.ini.h are incorrect')
        end
        if input.mdv.nf~=1
            if size(input.ini.Fak,1)~=input.mdv.nef
                error('The rows of input.ini.Fak need to be the number of size fractions minus 1')
            elseif size(input.ini.Fak,2)~=1 && size(input.ini.Fak,2)~=input.mdv.nx
                error('The columns of input.ini.Fak need to be the number of x points')
            end
        end
        %check slope or bed elevation
        if isfield(input.ini,'etab') && isfield(input.ini,'slopeb')
            error('If you provide the initial bed elevation you cannot provide the initial slope')
        end
        if (isfield(input.ini,'slopeb') && isfield(input.ini,'etab0')==0) || (isfield(input.ini,'slopeb')==0 && isfield(input.ini,'etab0'))
            error('If you provide the initial slope you have to provide the downstream bed elevation and viceversa')
        end
        if isfield(input.ini,'etab')
            if numel(input.ini.etab)~=1 && size(input.ini.etab,2)~=input.mdv.nx
                error('The dimensions of input.ini.etab are incorrect')
            end
        end
        if isfield(input.ini,'slopeb')
            if numel(input.ini.slopeb)~=1 && size(input.ini.slopeb,2)~=input.mdv.nx
                error('The dimensions of input.ini.slopeb are incorrect')
            end
        end

    case 3 %from file
        %check that inicon.mat exists
        warning('there is no check to the input when willing to load file')
        
    case 4 %normal flow (for a given initial condition)
        input.bch.mor='set1';
        input.bch.dotype='set1';
        input.bcm.type='set1';
        if isfield(input.ini,'u') || isfield(input.ini,'h')
            fprintf(fid_log,'ATTENTION!!! input.ini.u and input.ini.b are computed to have normal flow for a given water discharge and slope \n');
        end
        %check dimension if vectors
        if input.mdv.nf~=1
            if size(input.ini.Fak,1)~=input.mdv.nef
                error('The rows of input.ini.Fak need to be the number of size fractions minus 1')
            elseif size(input.ini.Fak,2)~=1 && size(input.ini.Fak,2)~=input.mdv.nx
                error('The columns of input.ini.Fak need to be the number of x points')
            end
        end
        
    case {5,51,52,53,12,13} %alternating steady equilibrium profile
        warningprint(fid_log,'For a space-marching initial condition no input-check is performed');
       
    otherwise
        error('input.ini.initype can be: 1 (normal flow), 2 (free), 3 (from file), 4 (normal flow out of initial condition), 5(space marching)')
end

   %check substrate
	if input.mdv.nf~=1
        if strcmp(input.ini.fsk,'Fak')==1
            warningprint(fid_log,'The substrate composition is assumed to be the same as the active layer composition');
        else
            %load it if it is a file
            if ischar(input.ini.fsk)
                load(input.ini.fsk);
            end
                %actual check
            if     size(input.ini.fsk,1)~=input.mdv.nef
                    error('The rows of input.ini.fsk need to be the number of size fractions minus 1')
            elseif size(input.ini.fsk,2)~=1 && size(input.ini.fsk,2)~=input.mdv.nx        
                    error('The columns of input.ini.fsk need to be the number of x points')
            elseif size(input.ini.fsk,3)~=1 && size(input.ini.fsk,3)~=input.mdv.nsl        
                    error('The columns of input.ini.fsk need to be the number of layers in the substrate')
            end
               %add patch
            if isfield(input.ini.fsk,'patch') 
                if any(input.ini.fsk.releta<0)
                    error('Specify substrate depth, positive downward')
                end
                %check dimension if vectors
                if size(input.ini.fsk,1)~=input.mdv.nef
                    error('The rows of input.ini.fsk need to be the number of size fractions minus 1')
                elseif size(input.ini.fsk,2)~=1
                    error('The columns of input.ini.fsk needs to be 1')
                end 
            end
        end
	end


%%
%% BCH
%%


%% UPSTREAM

%check that the kind of upstream hydrodynamic boundary condition is specified
if isfield(input.bch,'uptype')==0
    error('You need to specify the kind of upstream hydrodynamic boundary condition in input.bch.uptype')
end

%check that the input matches the kind of input
switch input.bch.uptype
    case {1,11} %water discharge Q
        if isfield(input.bch,'timeQ0')==0
            error('You need to provide the time at which the specific water discharge is specified (input.bch.timeQ0) if you want such an upstream hydrodynamic boundary condition')
        elseif isfield(input.bch,'Q0')==0
            error('You need to provide the specific water discharge at the specified times (input.bch.Q0) if you want such an upstream hydrodynamic boundary condition')
        end 
        if input.bch.timeQ0(end)<input.mdv.dt
            error('The last time of the upstream hydraulic boundary condition in the boundary condition is not larger or equal to the chosen time step')
        elseif input.bch.timeQ0(end)<input.mdv.Tstop
            warningprint(fid_log,'The last time of the upstream hydraulic boundary condition in the boundary condition is not larger or equal to the simulation time. The hydrograph will be repeated.')
        end

    case {12,13,14}
        if isfield(input.bch,'path_file_Q')==0 %if there does not exist a path to the Q-file 
            [path_folder_main,~,~]=fileparts(path_file_input); %get path to main folder
            input.bch.path_file_Q=fullfile(path_folder_main, 'Q.mat'); %path to the Q file
        end 

    case 'set1'
        if isfield(input.bch,'timeQ0') || isfield(input.bch,'Q0')
            fprintf(fid_log,'ATTENTION! input.bch.timeQ0 and input.bch.Q0 are not used if you specify normal flow for a given initial condition \n');
        end
        
        
    otherwise
        error('input.bch.uptype can be: 1, 11=water discharge from input, 12=water discharge from file')
end        



%% DOWNSTREAM

%check that the kind of downstream hydrodynamic boundary condition is specified
if isfield(input.bch,'dotype')==0
    error('You need to specify the kind of downstream hydrodynamic boundary condition in input.bch.dotype')
end

%check that the input matches the kind of input
switch input.bch.dotype
    case {1,11} %downstream water elevation
        if isfield(input.bch,'timeetaw0')==0
            error('You need to provide the time at which the downstream water level is specified (input.bch.timeetaw0) if you want such a downstream hydrodynamic boundary condition')
        elseif isfield(input.bch,'etaw0')==0
            error('You need to provide the donwstream water level at the specified times (input.bch.etaw0) if you want such an upstream hydrodynamic boundary condition')
        end 
        if input.bch.timeetaw0(end)<input.mdv.dt
            error('The last time of the downstream hydraulic boundary condition in the boundary condition is not larger or equal to the chosen time step')
        elseif input.bch.timeetaw0(end)<input.mdv.Tstop
            warningprint(fid_log,'The last time of the downstream hydraulic boundary condition in the boundary condition is not larger or equal to the simulation time. The downstream boundary condition will be repeated') 
        end
    case 12
        if isfield(input.bch,'path_file_etaw0')==0 %if there does not exist a path to the hdown-file 
            [path_folder_main,~,~]=fileparts(path_file_input); %get path to main folder
            input.bch.path_file_etaw0=fullfile(path_folder_main, 'etaw0.mat'); %path to the hdown file
        end         
    case 'set1'
        if isfield(input.bch,'timeetaw0') || isfield(input.bch,'etaw0')
            fprintf(fid_log,'ATTENTION!!! input.bch.timeetaw0 and input.bch.etaw are computed to have normal flow for a given water discharge and slope \n');
        end
    case {2,21} %downstream water depth
        error('not working yet??')
%         if isfield(input.bch,'path_file_etaw0')==0 %if there does not exist a path to the hdown-file 
%             [path_folder_main,~,~]=fileparts(path_file_input); %get path to main folder
%             input.bch.path_file_etaw0=fullfile(path_folder_main, 'etaw0.mat'); %path to the hdown file
%         end 
    case 22 %water depth from file
        error('not working yet??')
    otherwise
        error('input.bch.dotype can be: 1,11=downstream water elevation from input, 12=downstream water elevation from file')
end

%%
%% BCM
%%

%check that the kind of morphodynamic boundary condition is specified
if isfield(input.bcm,'type')==0
    error('You need to specify the kind of morphodynamic boundary condition in input.bcm.type')
end

%check that the input matches the kind of input
switch input.bcm.type
    case {1,11} %sediment discharge
        if isfield(input.bcm,'timeQbk0')==0
            error('You need to provide the time at which the volume of sediment transported excluding pores per unit time, and per size fraction is specified (input.bcm.timeQbk0) if you want such a morphodynamic boundary condition')
        elseif isfield(input.bcm,'Qbk0')==0
            error('You need to provide the volume of sediment transported excluding pores per unit time, and per size fraction at the specified times (input.bcm.Qbk0) if you want such a morphodynamic boundary condition')
        end 
        %check the dimensions of the input
        if size(input.bcm.Qbk0,1)~=size(input.bcm.timeQbk0,1)
            error('The rows of input.bcm.Qbk0 must be the same as the rows of input.bcm.timeQbk0')
        end
        if size(input.bcm.Qbk0,2)~=input.mdv.nf
            error('The colums of input.bcm.Qbk0 must be the same as the number grain size fractions')
        end
        %check proper time definition
        if input.bcm.timeQbk0(end)<input.mdv.dt
            error('The last time of the upstream morphodynamic boundary condition in the boundary condition is not larger or equal to the chosen time step')
        elseif input.bcm.timeQbk0(end)<input.mdv.Tstop
            warningprint(fid_log,'The last time of the upstream morphodynamic boundary condition in the boundary condition must be larger or equal to the simulation time')
        end


    case 12
        if isfield(input.bcm,'path_file_Qbk0')==0 %if there does not exist a path to the Qbk-file 
            [path_folder_main,~,~] = fileparts(path_file_input); %get path to main folder
            input.bcm.path_file_Qbk0=fullfile(path_folder_main, 'Qbk0.mat'); %path to the Qbk file
        end 

    case {2,'set1'}
	%nothing to check
    
    case 13
    %nothing implemented yet... 
    % check input.bcm.NFLtype;input.bcm.NFLparam; input.bcm.NFLiparam (at
    % least check the dimensions)
    
    
    otherwise
        error('input.bcm.type can be: 1=sediment discharge; 2=periodic; 12 ')
end

%%
%% MOR
%%

%active layer type
if input.mdv.nf~=1
    %interfacetype
    if input.mor.interfacetype==2
        if isfield(input.mor,'fIk_alpha')==0
            error('Specify fIk_alpha if you use Hoey and Ferguson')
        end
        if input.mor.fIk_alpha<0 || input.mor.fIk_alpha>1
            error('input.mor.fIk_alpha must be between 0 and 1')
        end
    end  
    switch input.mor.Latype
        case 4
            if isfield(input.mor,'La_t_growth')==0
                error('Specify the growing factor of the active layer')
            end
    end
else
    input.mor.Latype=1; %pseudo
end

%thickness of the last layer
if isfield(input.mor,'ThUnLyrEnd')==0 %if it is not defined 
    input.mor.ThUnLyrEnd=input.mor.ThUnLyr; %use the same as the rest
end

%morfac
if isfield(input.mor,'MorFac')==0 %if it does not exist
    input.mor.MorFac=1; 
end

%ellipticity check
if isfield(input.mor,'ellcheck')==0 %if it does not exist
    input.mor.ellcheck=0; %default is not to check (expensive!)
end

if input.mor.ellcheck==0 && input.mdv.chk.ell==1
    input.mor.ellcheck=1;
    warningprint(fid_log,'Hum... seems you want to check for ellipticity but you have not actually asked for it (check input.mor.ellcheck). I am checking for ellipticity because it seems you are a nice person, but ELV does not like stupid input! :D')
end

%% 
%% TRA
%% 

input=add_sedflags(input);

%%
%% NOUR
%%

if isfield(input,'nour')==0
    input.nour.t=NaN;
else
    if input.mdv.nf==1
        input.nour.dk_opt=NaN;
    end
end

%%
%% OUTPUT CHECK
%%

order_result=0;
if     any(strcmp(input.mdv.output_var,'u'))
    order_result=order_result+input.mdv.nx;
elseif any(strcmp(input.mdv.output_var,'h'))
    order_result=order_result+input.mdv.nx;
elseif any(strcmp(input.mdv.output_var,'etab'))
    order_result=order_result+input.mdv.nx;
elseif any(strcmp(input.mdv.output_var,'La'))
    order_result=order_result+input.mdv.nx;
elseif any(strcmp(input.mdv.output_var,'Cf'))
    order_result=order_result+input.mdv.nx;    
elseif any(strcmp(input.mdv.output_var,'ell_idx'))
    order_result=order_result+input.mdv.nx;    
elseif any(strcmp(input.mdv.output_var,'Mak'))
    order_result=order_result+input.mdv.nx*input.mdv.nef;    
elseif any(strcmp(input.mdv.output_var,'msk'))
    order_result=order_result+input.mdv.nx*input.mdv.nsl*input.mdv.nef;    
elseif any(strcmp(input.mdv.output_var,'Ls'))
    order_result=order_result+input.mdv.nx*input.mdv.nsl*input.mdv.nef;        
elseif any(strcmp(input.mdv.output_var,'qbk'))
    order_result=order_result+input.mdv.nx*input.mdv.nef;        
end
   
order_result=order_result*input.mdv.nT;

max_order_result=1e7;
warning_time=10; %[s]

warning('off','backtrace');
if order_result>max_order_result
    for ktw=warning_time:-1:0
        warning('Results file has %4.2e values, your computer may freeze. I will continue in %d seconds',order_result,ktw)
        pause(1)
    end
end
warning('on','backtrace');

%%
%% SAVE
%% 

save(fullfile(path_folder_main,'input.mat'),'input')
input_out=input;
