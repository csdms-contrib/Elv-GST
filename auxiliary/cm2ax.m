%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                       ELV                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%This awesome model has been created by Liselot and Victor.
%Please use it with a lot of care and love. If you have any
%problem send us an email:
%v.chavarriasborras@tudelft.nl
%
%$Revision: 49 $
%$Date: 2017-04-04 09:33:12 +0200 (Tue, 04 Apr 2017) $
%$Author: V $
%$Id: cm2ax.m 49 2017-04-04 07:33:12Z V $
%$HeadURL: https://131.180.60.193/svn/ELV/branches/V0123/auxiliary/cm2ax.m $
%
%transforms a position in an axes in centimeters to its value in axes units
%
%INPUT:
%   -pos_cm = desired position in cm (x,y); double [2,1]
%   -han_fig = handle to the figure; handle
%   -han_ax = handle to the axes; handle
%
%OUTPUT:
%   -pos_ax = position in the units of the axes (x,y); double [2,1]
%
%OPTIONAL:
%   -'reference' = reference point where the measure is given
%       -'ll' = lower left corner of the axes (Position)
%       -'lr' = lower right corner of the axes (Position)
%       -'ul' = upper left corner of the axes (Position)
%       -'ur' = upper right corner of the axes (Position)

function pos_ax=cm2ax(pos_cm,han_fig,han_ax,varargin)
%% parse
parin=inputParser;

input.reference.default='ll'; %lowerleft
addOptional(parin,'reference',input.reference.default);

parse(parin,varargin{:});

reference=parin.Results.reference;

%% axes dimensions

han_fig.Units='centimeters';
ax_width_cm=han_ax.Position(3)*han_fig.PaperPosition(3);
ax_height_cm=han_ax.Position(4)*han_fig.PaperPosition(4);

%% reference
switch reference
    case 'll'
        x_cm=pos_cm(1)/ax_width_cm;
        y_cm=pos_cm(2)/ax_height_cm;
    case 'lr'
        x_cm=(ax_width_cm-pos_cm(1))/ax_width_cm;
        y_cm=pos_cm(2)/ax_height_cm;
    case 'ul'
        x_cm=pos_cm(1)/ax_width_cm;
        y_cm=(ax_height_cm-pos_cm(2))/ax_height_cm;
    case 'ur'
        x_cm=(ax_width_cm-pos_cm(1))/ax_width_cm;
        y_cm=(ax_height_cm-pos_cm(2))/ax_height_cm;
    otherwise
        error('the reference can be: ll (lower left) | lr (lower right) | ul (upper left) | ur (upper right)')
end

%coordinates in axes dimensions
switch han_ax.XScale
    case 'linear'
        x_ax=x_cm*(han_ax.XLim(2)-han_ax.XLim(1))+han_ax.XLim(1);
    case 'log'
        x_ax=10^(x_cm*(log10(han_ax.XLim(2))-log10(han_ax.XLim(1)))+log10(han_ax.XLim(1)));
    otherwise
        error('The XScale can be: linear | log')
end
switch han_ax.YScale
    case 'linear'
        y_ax=y_cm*(han_ax.YLim(2)-han_ax.YLim(1))+han_ax.YLim(1);
    case 'log'
        y_ax=10^(y_cm*(log10(han_ax.YLim(2))-log10(han_ax.YLim(1)))+log10(han_ax.YLim(1)));
    otherwise
        error('The YScale can be: linear | log')
end

pos_ax=[x_ax,y_ax];
    
end