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
    
    model = @( p ) sum( (h_data - (offset + spc_model( p )).^2  ))/numel( data );
    %TODO: the debug arrived here: it seems that offset (mtf) and spc_model
    %      have different dimensions. Somehow.
    opts = optimset( 'MaxIter', 1e30, 'TolFun', 1e-30, 'TolX', 1e-30 );
    
    pees = fminsearch( model, pp, opts );
    while pp ~= pees
        pp = pees;
        pees = fminsearch( model, pp, opts );
    end
end
        
