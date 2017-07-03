function eta = elevation_reconstruction_mean(qb_down,qb_up,eta_mean,input,dx,dT)

% Total transport
qb_downt = sum(qb_down,3); 
qb_upt = sum(qb_up,3); 

% Construct time series bed elevation
dqbdx = -(qb_upt-qb_downt)/dx;
detadt = -1/(1-input.mor.porosity)*dqbdx*dT;

%eta = zeros(numel(qb_downt)+1,1);

%for j=1:numel(qb_downt);
%    dqbdx = -(qb_upt(j)-qb_downt(j))/dx;
%    detadt = -1/(1-input.mor.porosity)*dqbdx;
%    eta(j+1) = eta(j) + detadt*dT;
%end
eta = [0; cumsum(detadt)];
mean_eta = mean(eta);
eta = eta_mean + eta - mean_eta;

end