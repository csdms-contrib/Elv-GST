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
%$Date: 2017-04-20 13:07:54 +0200 (do, 20 apr 2017) $
%$Author: V $
%$Id: ini_spacem.m 79 2017-04-20 11:07:54Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/L0081/main/ini_spacem.m $
%
%active_reconstruction_mean ...
%
%\texttt{[eta, Fak] = active_reconstruction_mean(Qb_down,Qb_up,eta_mean,Fak_mean,input,dx,dT,B)}
%
%INPUT:
%   -
%
%OUTPUT:
%   -
%
%HISTORY:
%170613
%   -V. Forced Liselot to add a header

function [eta, Fak] = active_reconstruction_mean(Qb_down,Qb_up,eta_mean,Fak_mean,input,dx,dT,B,ic)   
    if isnan(Fak_mean)==1;
        Fak = NaN;
        eta_t0 = eta_mean;
        
        % Construct time series bed elevation
        dqbdx = -(Qb_up-Qb_down)/dx;        
        detadt = -1/(1-input.mor.porosity)*dqbdx*dT;
        
        % Get full series
        detadt = detadt(ic);        
        
        % Reconstruct
        %eta = [eta_t0; eta_t0+cumsum(detadt)];
        
        eta = [0; cumsum(detadt)];                   
        mean_eta = mean(eta);
        eta = eta_mean + eta - mean_eta;
        
        
        
    else
        % Total transport
        Qb_downt = sum(Qb_down,2); 
        Qb_upt = sum(Qb_up,2); 

        % Construct time series bed elevation
        dQbdx = -(Qb_upt-Qb_downt)/dx;
        detadt = -1/(B*(1-input.mor.porosity))*dQbdx*dT;

        % Get full series
        dQbdx = dQbdx(ic);
        detadt = detadt(ic);
        
        for j=1:numel(Fak_mean);
            fkI = Fak_mean(j);
            dQbkdx = -(Qb_up(:,j)-Qb_down(:,j))/dx;
            % Get full series
            dQbkdx = dQbkdx(ic);

            dFakdt = -fkI./input.mor.La(1,1).*detadt-1/(B*(1-input.mor.porosity))*dQbkdx*dT;
            Faktemp = [0; cumsum(dFakdt)];
            mean_Fak = mean(Faktemp);
            Fak(:,j) = Fak_mean(j) + Faktemp - mean_Fak;
        end
        eta = [0; cumsum(detadt)];                   
        mean_eta = mean(eta);
        eta = eta_mean + eta - mean_eta;

    end
end