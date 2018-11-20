%this utility calculates a numerical covariance matrix
%
% Cov = xb_covariance( function_handle, parameters )
%
% function_handle : a handle to a function tha takes a parameter vector
% parameters : the parameter vector to pass the function

function C = xb_covariance( fh, pees )
    C = inv( xb_hessian( fh, pees ) );
end
