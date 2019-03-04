%this utility redraws a spectrum from the output of xb_make_spc_ffb
%if the field array is given, a rebinning is possible.
%
% xb_redraw_ffb( hst_ffb, binz_ffb )
% [hst_ffb, binz_ffb, herr_ffb] = xb_redraw_ffb( nrg_ffb, bin, [figure] );
%
%parameters:
% --hst_ffb: the cell array containing the histograms
% --binz_ffb: the cell array containing the bins
% --nrg_ffb: the cell array with the content of the field
% --bin: a number or an array defining the bin
% --figure: optionally, a figure handle

function [hst_ffb, binz_ffb, herr_ffb] = xb_redraw_ffb( varargin )
    if nargout == 0 && iscell( varargin{1} ) && iscell( varargin{2} )
        hst_ffb = varargin{1};
        binz_ffb = varargin{2};
        fig = figure;
    elseif nargout > 0 && iscell( varargin{1} ) && nargin >= 2
        nrg_ffb = varargin{1};
        bin = varargin{2};
        if nargin == 3
            fig = varargin{3};
        end
    else
        error( 'Inconsistent arguments' );
    end
    
    if ~exist( 'hst_ffb' )
        hst_ffb = cell( 3, 1 ); binz_ffb = cell( 3, 1 );
        for ii=1:3
            [hst_ffb(ii), binz_ffb(ii)] = hist( nrg_ffb{ii}, bin );
        end
    end
            
    if exist( 'fig' ) && isfigure( fig );
		figure( fig );
		stairs( binz_ffb{1}, hst_ffb{1}, 'linewidth', 2 ); hold on;
		stairs( binz_ffb{2}, hst_ffb{2}, 'linewidth', 2 );
		stairs( binz_ffb{3}, hst_ffb{3}, 'linewidth', 2 );
		set( gca, 'fontsize', 24, 'linewidth', 2, 'yscale', 'log' );
		grid on;
		bstep = binz_ffb{1}(2)-binz_ffb{1}(1);
		ylabel( ['#/',num2str(bstep),'KeV'] );
		xlabel( 'KeV' );
		legend( { 'full CB', 'front CB', 'back CB' } );
	end
end
