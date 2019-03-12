%this utility takes in either a function handle or an histogram and provides a random
%number generator that returns numbers according to that distribution.
%
% arb_rng = xb_arbitrary_distro( function_handle, interval, [npart] );
% arb_rng = xb_arbitrary_distro( histogram, bin );
%
%arguments:
% A)
% -- function_handle : a function handle hopefully defined in interval.
% -- interval : a vector representing the interval. Can be also a binnage.
% -- npart : optionally, a partitioning for the interval (# slices).
% B)
% -- histogram : a collection of counts
% -- bin : the binnage of said histogram
%returns:
% -- arb_rng : a random generato (function handle) that provides numnbers
%              distributed according to the input.
%
%NOTE: The function handle, or histogram, is a pdf: it must be normalised on the interval
%      and it should not contain negative values. You'll be warned if these
%      conditions are not matched.
%      If the normalisation is _less_ than one no warning will be emitted, but you'll
%      get a significative number of NaNs. This is an inelegant way to convey information
%      though, so I'm not explicitly checking here.
%NOTE: These distros will give some very small percentage of NaNs: these are due to
%      the finite resolution of the sampling method (an histo with a finite amount of
%      bins or a function with a finite stepping). Generally, the smaller the bin the
%      less NaNs you'll get.

function arb_rng = xb_arbitrary_distro( distro, interval, npart )
    if is_function_handle( distro ) && nargin == 2
        arb_dist = __from_function( distro, interval );
    elseif is_function_handle( distro ) && nargin == 3
        if npart < 1
            error( 'npart must be at least 1' );
        end
        interval = linspace( interval(1), interval(2), npart );
        arb_dist = __from_function( distro, interval );
    elseif isvector( distro ) && nargin == 2
        arb_dist = __from_histo( distro, interval );
    else
        error( 'Inconsistent arguments: function handle and interval or histogram and binnage.' );
    end
    
    arb_rng = @( r, c ) arb_dist( rand( r, c, 'double' ) );
end

function arb_dist = __from_function( fh, interval )
    fsample = fh( interval );
    if ~isempty( find( fsample < 0 ) )
        warning( 'A pdf with negative values? Are YOU sure?' );
    end
    arb_dist = __from_histo( fsample, interval );
end

function arb_dist = __from_histo( h, b )
    if sum( h ) > 1
        warning( 'Your pdf is not normalised: YOU will not be able to access part of it.' );
    end
    m = tril( ones( numel( h ) ) );
    hs = m*h(:);
    
    arb_dist = @(flat) interp1( hs, b, flat );
end
