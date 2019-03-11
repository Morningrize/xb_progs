% An utility to rescale spectra from counts to reaction probability.
%
% [h_rp, herr_rp] = xb_counts2rp( scalers, h )
% [h_rp, herr_rp] = xb_counts2rp( scalers, h, b, fig )
%
%arguments:
% --scalers: the scalers structure returned by xb_tpat2scaler
% --h: an histogram or array of histograms
% --b: optionally, the binnage (to plot)
% --figure: optionally, a figure where to plot the historgram
%returns:
% --h_rp: the scaled spectrum
% --herr_rp: the poisson-uncertainties on the bins.
%
%NOTE: this makes sense only on already cut data. Before that, saying thing about
%      scaling --> reaction probability for s412 is as good as impossible.

function [h_rp, herr_rp] = xb_counts2rp( scalers, h, varargin )
    if iscell( h )
        h_rp = cell( size( h ) );
        herr_rp = cell( size( h ) );
        for ii=1:length( h )
            h_rp(ii) = h{ii}*scalers.a/sum( h{ii} );
            herr_rp = sqrt( h_rp{ii} );
        end
    elseif isvector( h )
        h_rp = h*scalers.a/sum( h );
        herr_rp = sqrt( h_rp );
    else
        error( 'this histogram does not make sense.' );
    end

    if nargin == 4 && isfigure( varargin{2} )
        fig = varargin{2};
        b = varargin{1};
        figure( fig );
        
        if iscell( h_rp )
            xb_redraw_ffb( h_rp, b, fig );
            bstep = b{1}(2) - b{1}(1);
        else
            stairs( b, h_rp, 'linewidth', 2 );
            bstep = b(2) - b(1);
        end
        ylabel( ['rp/',num2str(bstep),'KeV'] );
    end
end
