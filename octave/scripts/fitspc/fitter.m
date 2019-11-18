%This function does the fitting of a spectrum on top of a dataset.
%Usually the spectrum comes out of simulation and the data from a detector.
%
% [pees, pee_errs, chisq] = fitter( spc_pees, spc_model, h_data, extremes, binZ, minopts )
%
%parameters:
% spc_pees: the parameters to build the hybrid spectrum
% spc_model: a function handle that takes the parameters and returns a spectrum
% h_data : the spectrum (already histogram!) of the data.
% extremes : the interval where to perform the fit.
%returns
% pees : the result of the fit.
% pee_err : the errors thereupon (still experimental)
% chisq : the reduced chi squared of the fit

function [pees, pee_errs, chisq] = fitter( spc_pees, spc_model, h_data, extremes, ...
                                           binZ, minopt )
    L_h_data = log( max( h_data, 1 ) );
    L_spc_model = @( p ) log( max( spc_model( p ), 1 ) );
    
    if ~exist( 'extremes' )
        model = @( p ) sum( norm(L_h_data - L_spc_model( p )) )/numel( h_data );
    else
        model = @( p ) sum( norm( L_h_data(extremes(1):extremes(2)) - ...
                                  L_spc_model( p )(extremes(1):extremes(2)) ) )/ ...
                                  numel( extremes(1):extremes(2) );
    end
    
    pcage = [ones( 1, numel( spc_pees ) );zeros( 1, numel( spc_pees ) )];
    if exist( 'minopt', 'var' )
        [pees, jval, rc] = xb_constrained_gradient_descent( model, spc_pees, pcage, minopt );
    else
        [pees, jval, rc] = xb_constrained_gradient_descent( model, spc_pees, pcage );
    end
    
    %NOTE: this is still experimental and might very well be BS
    J_cov = xb_covariance( model, pees );
    pee_errs = abs( sqrt( diag( J_cov )' )/sqrt( extremes(2) - extremes(1) ) );
    chisq = xb_goodness_of_fit( h_data(extremes(1):extremes(2)), ...
                                spc_model( pees )(extremes(1):extremes(2)),...
                                numel( pees ) );
end
        
