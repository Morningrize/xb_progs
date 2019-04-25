% An utility to rescale spectra from reaction probability to counts.
%
% [h_cnt, herr_cnt] = xb_counts2rp( scalers, h, total )
% [h_cnt, herr_cnt] = xb_counts2rp( scalers, h, total, b, fig )
%
%arguments:
% --scalers: the scalers structure returned by xb_tpat2scaler
% --h: an histogram or array of histograms
% --b: optionally, the binnage (to plot)
% --figure: optionally, a figure where to plot the historgram
%returns:
% --h_cnt: the scaled spectrum
% --herr_cnt: the poisson-uncertainties on the bins.
%
%NOTE: this makes sense only on already cut data. Before that, saying thing about
%      scaling --> reaction probability for s412 is as good as impossible.

function [h_cnt, herr_cnt] = xb_rp2counts( scalers, h, total, varargin )
    if iscell( h )
        h_cnt = cell( size( h ) );
        herr_cnt = cell( size( h ) );
        for ii=1:length( h )
            h_cnt(ii) = h{ii}/scalers.a*total;
            herr_cnt = sqrt( h_cnt{ii} );
        end
    elseif isvector( h )
        h_cnt = h/scalers.a*total;
        herr_cnt = sqrt( h_cnt );
    else
        error( 'this histogram does not make sense.' );
    end

    if nargin == 4 && isfigure( varargin{2} )
        fig = varargin{2};
        b = varargin{1};
        figure( fig );
        
        if iscell( h_cnt )
            xb_redraw_ffb( h_cnt, b, fig );
            bstep = b{1}(2) - b{1}(1);
        else
            stairs( b, h_cnt, 'linewidth', 2 );
            bstep = b(2) - b(1);
        end
    end
end
