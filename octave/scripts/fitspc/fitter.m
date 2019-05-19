%This function does the fitting of a spectrum on top of a dataset.
%Usually the spectrum comes out of simulation and the data from a detector.
%
% fitted_parameters = fitter( spc_pees, spc_model, h_data, extremes, binZ, offset )
%
%parameters:
% spc_pees: the parameters to build the hybrid spectrum
% spc_model: a function handle that takes the parameters and returns a spectrum
% h_data : the spectrum (already histogram!) of the data.
% extremes : the interval where to perform the fit.
% offset : an eventual offset spectrum
%returns
% fitted_parameters : the result of the fit.

function pees = fitter( spc_pees, spc_model, h_data, extremes, offset )
    pp = spc_pees;
    pees = zeros( size( spc_pees ) );

    model = @( p ) sum( (h_data - spc_model( p ) - offset).^2 )/numel( h_data );
    opts = optimset( 'MaxIter', 1e9, 'TolFun', 1e-9, 'TolX', 1e-9 );
    
    ii=0;
    pees = fminsearch( model, pp, opts );
    while pp ~= pees
        pp = pees;
        pees = fminsearch( model, pp, opts );

        screwed = find( pees > 1e3 );
        if ~isempty( screwed )
            warning( ['There was a runaway at iteration',num2str(ii)] )
            pees( screwed ) = rand( numel( screwed ) );
        end
        
        ++ii;
        if ii > 1e3
            warning( 'Could not converge in 1000 attempts' );
            break;
        end
    end
    
    %uncertainties estimation missing still!
end
        
