%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 120 $
%$Date: 2017-06-29 10:13:24 +0200 (Thu, 29 Jun 2017) $
%$Author: V $
%$Id: flow_update.m 120 2017-06-29 08:13:24Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/main/flow_update.m $
%
%flow_update updates the flow depth and the mean flow velocity.
%
%[u_new,h_new]=flow_update(u,h,etab,Cf,bc,input,fid_log,kt)
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
%160418
%   -L. bwc bug with qw=q_old
%

%160428
%   -L. Update hbc and Qbc cases
%       - Update case numbers
%       - Modulo construction for repeated time steps
%

%160513
%   -L. Implicit Preismann scheme
%
%
%160623
%   -V. cyclic boundary conditions
%
%160702
%   -V. tiny performance improvements
%
%160803
%	-L. Merged Vv4 Lv4
%
%160818
%   -L. Adjusted for variable flow.
%
%161130
%   -L. Geklooi met breedte_variaties

function [u_new,h_new]=flow_update(u,h,etab,Cf,bc,input,fid_log,kt)

%%
%% RENAME
%%


H_old=h'; %water depth [m]; [nxx1 double]
q_old=(u.*h)'; %specific water discharge [m^2/s]; [nxx1 double]

eta_old=etab'; %!!! do not change! some commented parts have it

nx=input.mdv.nx; %number of cells

flowtype=input.mdv.flowtype;

g=input.mdv.g;
dx=input.grd.dx;
dt=input.mdv.dt;

B = input.grd.B;

Qbct = mod(kt,bc.repQT(2))+(mod(kt,bc.repQT(2))==0)*bc.repQT(2);
hbct = mod(kt,bc.rephT(2))+(mod(kt,bc.rephT(2))==0)*bc.rephT(2);

%%
%% BOUNDARY CONDITION
%%

%% DOWNSTREAM

switch input.bch.dotype
%%water elevation    
    case {1,11,12,13,'set1'} 
%         Hdown=bc.etaw0(hbct)-etab(K); %water depth at the downsteam end [1x1 double] (ERROR! etab is at cell center and boundary condition at cell edge)
        Hdown=bc.etaw0(hbct)-3/2*etab(nx)+1/2*etab(nx-1); %water depth at the downsteam end [1x1 double]
%%water depth   
    case {2,21,22,23} 
        error('not yet implemented')
    otherwise
		error('The magic downstream boundary condition is not yet implemented. Please provide a input.bch.dotype that is implemented.')

end

%% UPSTREAM 

switch input.bch.uptype
%%water discharge    
    case {1,11,12,13,14,2}
        qwup=bc.q0(Qbct); %[1x1 double]
        Qwup=bc.Q0(Qbct); %[1x1 double]
	otherwise
        error('The answer is 42, why do yo need to run a model? In any case, please provide a input.bch.uptype that is implemented.')
end


%%
%% UPDATE
%%
switch flowtype
    
%% CONSTANT FLOW
    case 0
        U=u;
        H=h;
        
%% STEADY FLOW

    case 1
        %downstream condition
        H(nx,1)=Hdown;
        qw=(qwup*B(1)./B)';
        Energy(nx,1)=eta_old(nx,1)+Hdown+0.5*qw(nx,1)^2/(g*Hdown^2);
        Sf(nx,1)=Cf(1,nx)*(qw(nx,1)/(H(nx,1)))^2*(1/(g*H(nx,1))); %((1/sqrt(Cf_q(j,1)))*H(j,1)/qw(j,1)).^(-2)*1/(g*H(j,1)); 
        for kx=nx-1:-1:1           
            %  Compute the dimensionless friction term, Euler explicit towards the upstream direction
            Energy(kx,1)=Energy(kx+1,1)+dx*Sf(kx+1,1);
            %  Invert the energy function analytically
            acoeff=1;
            bcoeff=-(Energy(kx,1)-eta_old(kx,1));
            ccoeff=0;
            dcoeff=0.5*qw(kx,1)^2/g;
            root=analytical_cubic_root(acoeff,bcoeff,ccoeff,dcoeff);
            H(kx,1)=max(root); %Select the subcritical value
            Sf(kx,1)=Cf(1,kx)*(qw(kx,1)/(H(kx,1)))^2*(1/(g*H(kx,1))); %((1/sqrt(Cf_q(j,1)))*H(j,1)/qw(j,1)).^(-2)*1/(g*H(j,1)); 
        end
    U=qw./H; %[nxx1 double]
        
%% QUASI STEADY FLOW 
    
    case 2

        % Preismann scheme;
        % Implicit computation
        
        [U,H] = preissmann(u,h,etab,Cf,Hdown,qwup,input,fid_log,kt);
           
%% UNSTEADY FLOW

    case 3       

        %----------------- GHOST CELLS AROUND BOUNDARIES -----------------%
        
        % Store old data in Kx2 array
        Q = [H_old,q_old]; %size K
      
        % Adjust BCs
        Q(1,2) = qwup;
        Q(end,1)= Hdown;
        
        % Q has 4 ghost cells; two on each end.       
        % First order of extrapolation at both ends
        Qgu1 = Q(1,:) - (Q(2,:)-Q(1,:));
        Qgu2 = Qgu1 - (Q(1,:)-Qgu1);
        Qgd1 = Q(end,:) + (Q(end,:)-Q(end-1,:));
        Qgd2 = Qgd1 + (Qgd1-Q(end,:));
        Q = [Qgu2; Qgu1; Q; Qgd1; Qgd2];
       
        % -------------- COMPUTATION AQ+ and AQ- ------------------------%
        % The system will be solved by summing the contribution of the
        % individual wave is in the system. Therefore the eigenvalues and
        % left and right eigenvectors are computed. From now on, 1 refers
        % to the positive eigenvalue, and 2 to the negative one.
        
        l11= Q(:,2)./Q(:,1) + sqrt(g*Q(:,1));
        l22= Q(:,2)./Q(:,1) - sqrt(g*Q(:,1));
        
        temp = 2*sqrt(g*Q(:,1));
        left_eig1 = [-l11.*l22./temp, l11./temp];
        left_eig2 = [l11.*l22./temp, -l22./temp];       
        right_eig1(:,:) = [1./l11(:,1)'; ones(1,nx+4)];
        right_eig2(:,:) = [1./l22(:,1)'; ones(1,nx+4)];
             
        Dq = Q(2:end,:)-Q(1:end-1,:); %difference between adjacent cells
        
        % Computation alpha; alpha1 is for the first eigenvalue (+), and
        % alpha2 for the second eigenvalue (-). Alpha is defined at cell
        % boundaries as the left eigenvalue multiplied by the difference in
        % Q. The upwind value for the left eigenvalue is used.
        
        
        % OLD LOOP :)
        %alpha1 = zeros(1,K+3); alpha2 = zeros(1,K+3);                       
        %parfor i=1:K+3
        %    alpha1(i) = left_eig1(i,:)*Dq(i,:)'; 
        %    alpha2(i) = left_eig2(i+1,:)*Dq(i,:)';  
        %end
        
        alpha1 = sum(left_eig1(1:nx+3,:).*Dq(1:nx+3,:),2);
        alpha2 = sum(left_eig2(2:nx+4,:).*Dq(1:nx+3,:),2);
        alpha1 = alpha1';
        alpha2 = alpha2';
                
        % OLD LOOP :)
        %AplusDq = zeros(K+3,2); AminDq = zeros(K+3,2);
%!!!!!! % For now compute in loop, replace later on;
        %for i=1:K+3
        %   AplusDq(i,:) = l11(i)*alpha1(i)*right_eig1(:,i); 
        %   AminDq(i,:) = l22(i+1)*alpha2(i)*right_eig2(:,i+1); 
        %end
        
        AplusDq = (ones(2,1)*(l11(1:nx+3)'.*alpha1(1:nx+3))).*right_eig1(:,1:nx+3);
        AminDq =  (ones(2,1)*(l22(2:nx+4)'.*alpha2(1:nx+3))).*right_eig2(:,2:nx+4);
                
        AplusDq = AplusDq';
        AminDq = AminDq';
        
        
        % ------------------ FLUX LIMITER CORRECTION STEP-----------------%
        % Computation of theta, theta at a cell boundary is defined as the
        % ration between alpha at the upwind cell boundary, divided by
        % alpha at the current cell boundary.       
        % Theta is size K+2
        theta1 = alpha1(1:end-1)./alpha1(2:end); %theta1(i-1/2)=alpha1(i-3/2)/alpha1(i-1/2)
        theta2 = alpha2(2:end)./alpha2(1:end-1); %theta2(i-1/2)=alpha2(i+1/2)/alpha2(i-1/2)
                
        % Compute flux limited alpha;
        alpha1_tilde = phi_func(theta1(1:nx+1),input).*alpha1(2:nx+2); %size K+1 
        alpha2_tilde = phi_func(theta2(2:nx+2),input).*alpha2(2:nx+2); %size K+1
        
        % OLD LOOP:)
        %Fl =zeros(K+1,2);
        %for i=1:K+1
        %    Fl(i,:) = 0.5*abs(l11(i+1)).*(1-dt/dx*abs(l11(i+1)))*alpha1_tilde(i)*right_eig1(:,i+1) + ...
        %    0.5*abs(l22(i+2)).*(1-dt/dx*abs(l22(i+2)))*alpha2_tilde(i)*right_eig2(:,i+2);              
        %end
        Fl = (ones(2,1)*(0.5*abs(l11(2:nx+2)').*(ones(1,nx+1)-dt/dx*abs(l11(2:nx+2)').*alpha1_tilde(1:nx+1)))).*right_eig1(:,2:nx+2) + ...
             (ones(2,1)*(0.5*abs(l22(3:nx+3)').*(ones(1,nx+1)-dt/dx*abs(l22(3:nx+3)').*alpha2_tilde(1:nx+1)))).*right_eig2(:,3:nx+3);
                
        
        
        % ----------------- SOURCE TERMS ---------------------------------%
        eta_old = [2*eta_old(1,1)-eta_old(2,1); eta_old; 2*eta_old(end,1)-eta_old(end-1,1)];
        S0 = -(eta_old(3:end,:) - eta_old(1:end-2,:))/(2*dx);
        S = [zeros(nx,1), g*Q(3:end-2,1).*S0 - (Cf' .* (Q(3:end-2,2).^2)) ./ (Q(3:end-2,1).^2)];
          
        % ----------------- UPDATE ---------------------------------------%
        Q_new = Q(3:nx+2,:) -(dt/dx)*(AplusDq(2:nx+1,:)+AminDq(3:nx+2,:)) +dt*S -(dt/dx)*(Fl(:,2:nx+1)'-Fl(:,1:nx)');

        qw = Q_new(:,2);
        H = Q_new(:,1);
        U = qw./H;


    case 4
        % Preismann scheme;
        % Implicit computation
        
        [U,H] = preissmann(u,h,etab,Cf,Hdown,qwup,input,fid_log,kt);
    
    case 5
        % already taken for space marching
        
    case 6
        % Backwater solver; 
        if input.grd.crt==1
            bedslope = [NaN; -(eta_old(2:end)-eta_old(1:end-1))/input.grd.dx];
            [U,H] = backwater_rect(bedslope,Cf,Hdown,Qwup,input);
        else
            bedslope = [NaN; -(eta_old(2:end)-eta_old(1:end-1))/input.grd.dx];
            [U,H] = backwater(bedslope,Cf,Hdown,Qwup,input);
        end
        
    otherwise
        error('We have not yet implemented the super flow solver that reads your mind to know how to solve the flow. Please provide a flowtype solver.')
end

%%
%% RENAME
%%
        u_new=U';
        h_new=H';
        


end

 
 

