%So, this is going to be a slow but hopefully reliable minimizer based on the good ol' gradient
%descent. This is the contrained version, bovine level of intelligence.
%
% [local_min, Jval] = xb_gradient_descent( function_handle, parameters, constraint, [options] )
%
%parameters:
% function_handle: a function handle that takes as parameters "parameters", and needs
%                  to be minimized.
% parameters: the initial guess for the parameters.
% constraint: the constraint to obey.
%             2*length( parameters ): first row --> upper; second row -- lower.
%retrurns:
% local_min: the local minimum (gradient descent only does that)
% Jval: the value of the function there.
%NOTE: there are little safery nets in here. Screw up the learning rate or the zero AND
%      you'll spin around indefinitely.

function [local_min, Jval, rc] = ...
    xb_constrained_gradient_descent( function_handle, parameters, constraint, varargin )
    lr = 1e-2;
    zr = 1e-9;
    mi = 1e5;
    rc = 0;
    
    for ii=1:2:numel( varargin )
        switch( varargin{ii} )
            case { 'lr', 'learning-rate' }
                lr = varargin{ ii+1 };
            case { 'z', 'zero-is' }
                zr = varargin{ ii+1 };
            case { 'M', 'max-iter' }
                mi = varargin{ ii+1 };
            otherwise
                warning( ['Unkown option: ',varargin{ii}] );
        end
    end
    
    ii=0;
    do
        dp = parameters;
        grad = xb_gradient( function_handle, parameters );
        parameters -= lr*grad;
        iup = find( parameters > constraint(1,:) );
        idn = find( parameters < constraint(2,:) );
        parameters( iup ) = constraint(1,iup);
        parameters( idn ) = constraint(2,idn);
        dp -= parameters;
        if( ii > mi )
            rc = 1;
            break;
        end;
        ++ii;
    until norm( grad ) < zr || abs( dp ) < zr
    
    local_min = parameters;
    Jval = function_handle( local_min );
end
