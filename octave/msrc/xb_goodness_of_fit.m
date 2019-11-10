%This function calculated the reduced chi squared of a fit
%
% chisq = xb_goodness_of_fit( data, model, [npees] )
%
%parameters:
% -- data: it's the data (histogram) that has been fitted
% -- model: it's the fitted model
% -- npees: optionally, the number of fit parameters
%returns
% -- chisq: the chi squared.

function chisq = xb_goodness_of_fit( data, model, npees )
    if numel( data ) ~= numel( model )
        error( 'Data and model must have the same numel' );
    end
    iv = find( model );
    chisq = sum( ( data(iv) - model(iv) ).^2./model(iv) );
    if ~exist( 'npees', 'var' ); npees = 0; end
    chisq /= numel( iv ) - npees;
end
