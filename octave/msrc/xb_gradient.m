%athis utility calculates the numerical gradient fo a function (handle) R^n -> R^m
%
% gradient = xb_gradient( function_handle, X0, delta );
%
%arguments:
% -- function_handle : is a function handle that accepts a vector as an argument.
% -- X0 : is the point at which the gradient will be calculated.
% -- delta : optionally, it's how small the step will be. Default is 1e-5.

function grad = xb_gradient( fh, x0, delta )
    if ~exist( 'delta', 'var' ); delta = 1e-5; end

    %dx can also be calculated relatively, but it should not be necessary!
    %dx = (x0*delta)(:);
    dx = delta*ones( 1, length( x0 ) );

    bunch = diag( x0 )*ones( length( x0 ), length( x0 ) );
    bunchplus = bunch + diag( dx );
    bunchminus = bunch - diag( dx );
    
    grad = [];
    for ii=1:length( x0 )
        grad = [grad, ( double( fh( bunchplus(:,ii) ) ) - double( fh( bunchminus(:,ii) ) ) )(:)/(2*dx(ii))];
    end
end
