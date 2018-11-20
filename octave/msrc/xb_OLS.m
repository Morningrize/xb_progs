%This function calculates the Ordinari Least Square regression (polynomial, interpreted as linear)
%for a given dataset and to a specified order. It also provides an estimate of the errors on the coefficients.
%
% [pees, pees_err] = xb_OLS( Xs, Ys, degree )
%
% Parameters:
% -- Xs: the independent variable.
% -- Ys: the measured values.
% -- degree: optionally, the degree of the polynomial to use --if omitted, defaults to 1
% Returns:
% -- pees: the fitted parameters
% -- pees_err: the errors thereupon.

function [ pees, pees_err ] = xb_OLS( Xs, Ys, degree )
    if nargin == 2
        degree = 1;
    end
    
    Xmat = ones( size( Xs(:) ) );
    
    for ii=1:degree
        Xmat = [Xmat, Xs(:).^ii];
    end
    
    pees = inv( Xmat'*Xmat )*Xmat'*Ys(:);
    
    if nargout == 1 return; end
    
    model = @( p ) sqrt( sum( (Ys - __poly( p, Xs )).^2 ) );
    pees_err = sqrt( diag( xb_covariance( model, pees ) ) );
end

function val = __poly( coeff, x )
    x = x(:);
    val = coeff(1)*ones( size(x) );
    for ii=2:length(coeff)
        val += coeff(ii)*x.^(ii-1);
    end
end
