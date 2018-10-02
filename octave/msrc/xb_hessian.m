%this utility calculate the hessian of a function
%
% H = xb_hessian( function_handle, x0 )
%
% function_handle : it's a function function handle
% x0: it's the point at which the hessian will be calculated

function H = xb_hessian( fh, x0 )
    gh = @( p ) xb_gradient( fh, p );
    H = xb_gradient( gh, x0 );
end
    

