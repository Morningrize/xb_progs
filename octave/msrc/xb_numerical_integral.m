%This utility performs a Monte Carlo numerical integration on a given function on a given interval
%NOTE: the function must be readily evaluable on the interval to be passed as a handle, otherwise
%      you're gonna need the evaluation.
%
% integral = xb_numerical_intgral( fh, interval, [nthrows] )
%
%parameters:
% -- fh: a function handle OR the function evaluation over the interval
% -- interval: EITHER the extremes of the interval
%              OR all the points where the function has been evaluated.
% -- nthrows: the number of dice throws the MC will use. defaults to 1e7,
%             which will give you a ~2e-4 precision on the 0.
%returns:
% -- integral: the value of the integral, as the MC generator think it should be.
%
%NOTE: the MC method used here is essentially based on a rectangle defined by the interval and
%      the function max value and then shooting at it very many times to get a probability of hit
%      that, multiplied by the area, will give the integral

function integral = xb_numerical_integral( fh, interval, nthrows )
    if ~exist( 'nthrows', 'var' ); nthrows = 1e7; end
    if is_function_handle( fh ) && numel( interval ) > 2
        interspan = interval;
        fh_eval = fh( interspan );
        fh = @(p) interp1( interspan, fh_eval, p );
        interval = [interval(1), interval(end)];
    elseif ~is_function_handle( fh ) && numel( interval ) == 2
        interspan = linspace( interval(1), interval(2), numel( fh ) );
        fh_eval = fh;
        fh = @(p) interp1( interspan, fh_eval, p );
    elseif ~is_function_handle( fh ) && numel( interval ) > 2
        if numel( fh ) ~= numel( interval ); error( 'fh has not been sampled on the interval' );
        else
            fh_eval = fh;
            interspan = interval; interval = [interval(1), interval(end)];
            fh = @(p) interp1( interspan, fh_eval, p );
        end
    end
    if ~exist( 'fh_eval', 'var' ); fh_eval = fh( linspace( interval(1), interval(2), 1e5 ) ); end

    idx_pos = find( fh_eval >= 0 );
    idx_neg = find( fh_eval < 0 );
    fh_pos = fh_eval( idx_pos );
    fh_neg = fh_eval( idx_neg );

    throws = rand( 2, nthrows );

    base = interval(2) - interval(1);
    height_pos = max( fh_pos ) - min( fh_pos );
    x_sampling = throws(1,:)*base + interval(1);
    y_sampling = throws(2,:)*height_pos + min( fh_pos );
    nb_pos = numel( find( y_sampling <= max( 0, fh( x_sampling ) ) ) );

    if isempty( fh_neg ); integral = nb_pos*nthrows^-1; return; end
    fh_neg *= -1;
    height_neg = max( fh_neg ) - min( fh_neg );
    y_sampling = throws(2,:)*height_neg + min( fh_neg );
    nb_neg = numel( find( y_sampling <= max( 0, -1*fh( x_sampling ) ) ) );

    integral = (nb_pos-nb_neg)*nthrows^-1;
end
