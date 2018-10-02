%this utility calculates the numerical gradient fo a function (handle) R^n -> R^m
%
% gradient = xb_gradient( function_handle, X0 );
%
%arguments:
% - function_handle is a function handle that accepts a vector as an argument.
% - X0 is the point at which the gradient will be calculated.

function grad = xb_gradient( fh, x0 )
    bunch = diag( x0 )*ones( length( x0 ), length( x0 ) );
    bunchplus = bunch + eye( length( x0 ) )*1e-5;
    bunchminus = bunch - eye( length( x0 ) )*1e-5;
    
    grad = [];
    for ii=1:length( x0 )
        grad = [grad, ( fh( bunchplus(:,ii) ) - fh( bunchminus(:,ii) ) )(:)/2e-5];
    end
end
