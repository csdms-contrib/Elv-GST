function obj = solve_nfbc(X,input,Qw,AL)
%solve_sedigraph computes the slope and surface fraction 
% VERSION 1

%INPUT:
%   -X(1) for slope, remainder for fractions
%   -input for parmaters
%   -Qw: equidistant spaced hydrograph
%   -AL: mean annual load per fraction
%
%OUTPUT:
%   -

%HISTORY:
%
%1600901
%   L-First creation

% Compute normal flow depth
Qw = (Qw/input.grd.B)';
ib = X(1);
h = (input.mdv.Cf.*Qw.^2/(9.81*ib)).^(1/3);

% Get lengths
nm = length(h);
nf = length(X);

% Initialize fractions
if nf>1
	%F1 = 1-sum(X(2:end));
	%Mak = repmat([F1;X(2:end-1)'],1,nm);
	Mak = repmat(X(2:end)',1,nm);
	La = ones(1,nm);
	if min(Mak)<0
    		warning('Negative fractions are not allowed. Revisit solver');
	end

else
	Mak = NaN*ones(size(h));
	La = NaN*ones(size(h));
end

Cf=input.mdv.Cf.*ones(size(h));  
if isfield(input,'tra.calib')==1
    [qbk,~]=sediment_transport(input.aux.flg,input.aux.cnt,h,Qw,Cf,La,Mak',input.sed.dk,input.tra.param,input.aux.flg.hiding_parameter,input.tra.calib,NaN,NaN);
else
    [qbk,~]=sediment_transport(input.aux.flg,input.aux.cnt,h,Qw,Cf,La,Mak',input.sed.dk,input.tra.param,input.aux.flg.hiding_parameter,1,NaN,NaN);
end

             
% Objectives:
% - Total mass per fraction is transported
% - Sum of all fractions is one. 
%pQ = ones(size(Qw))./numel(Qw);
%Computed_load = sum(qbk.*repmat(pQ,1,nf));
Computed_load = mean(qbk,1);
Annual_load = AL/input.grd.B;
if numel(AL) == 1
    obj = sum(Computed_load)-sum(Annual_load);
else    
    obj = Computed_load-Annual_load;
end
end

