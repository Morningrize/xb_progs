%This function does the fitting of a spectrum on top of a dataset.
%Usually the spectrum comes out of simulation and the data from a detector.
%
% [pees, pee_errs] = fitter( spc_pees, spc_model, h_data, extremes, binZ, offset )
%
%parameters:
% spc_pees: the parameters to build the hybrid spectrum
% spc_model: a function handle that takes the parameters and returns a spectrum
% h_data : the spectrum (already histogram!) of the data.
% extremes : the interval where to perform the fit.
% offset : an eventual offset spectrum
%returns
% pees : the result of the fit.
% pee_err : the errors thereupon (still experimental)

function [pees, pee_errs] = fitter( spc_pees, spc_model, h_data, extremes, offset )
    pp = spc_pees;
    pees = zeros( size( spc_pees ) );

    model = @( p ) sum( (h_data - spc_model( p ) - offset).^2 )/numel( h_data );
    opts = optimset( 'MaxIter', 1e12, 'TolFun', 1e-12, 'TolX', 1e-12 );
    
    ii=0;
    pees = fminunc( model, pp, opts );
    while sum( (pp - pees).^2 ) > 1e-15
        pp = pees;
        pees = fminunc( model, pp, opts );

        screwed = find( pees < 0 );
        if ~isempty( screwed )
            warning( ['There was a runaway at iteration ',num2str(ii)] )
            pees( screwed ) = rand( numel( screwed ), 1 );
        end
        
        ++ii;
#         if ii > 1e3
#             warning( 'Could not converge in 1000 attempts' );
#             break;
#         end
    end
    ii
    
    %NOTE: this is still experimental and might very well be BS
    J_cov = xb_covariance( model, pees );
    pee_errs = sqrt( diag( J_cov ) );
end
        
