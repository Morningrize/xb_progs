%So, this is going to be a slow but hopefully reliable minimizer based on the good ol' gradient
%descent.
%
% [local_min, Jval, rc] = xb_gradient_descent( function_handle, parameters, [options] )
%
%parameters:
% function_handle: a function handle that takes as parameters "parameters", and needs
%                  to be minimized.
% parameters: the initial guess for the parameters.
%retrurns:
% local_min: the local minimum (gradient descent only does that)
% Jval: the value of the function there.
% rc : the exit status:
%      -- 0 : reached convergence
%      -- 1 : hit the max iteration limit
%NOTE: there are little safery nets in here. Screw up the learning rate or the zero AND
%      you'll spin around indefinitely.

function [local_min, Jval, rc] = xb_gradient_descent( function_handle, parameters, minopts )
    lr = 1e-2;
    zr = 1e-9;
    mi = 1e5;
    rc = 0;
    
    if nargin == 3
        for ii=1:2:numel( minopts )
            switch( minopts{ii} )
                case { 'lr', 'learning-rate' }
                    lr = minopts{ ii+1 };
                case { 'z', 'zero-is' }
                    zr = minopts{ ii+1 };
                case { 'M', 'max-iter' }
                    mi = minopts{ ii+1 };
                otherwise
                    warning( ['Unkown option: ',minopts{ii}] );
            end
        end
    end
    
    ii = 0;
    do
        dp = parameters;
        grad = xb_gradient( function_handle, parameters );
        parameters -= lr*grad;
        if( ii > mi )
            rc = 1;
            break;
        end;
        ++ii;
        dp -= parameters;
    until norm( grad ) < zr || abs( dp ) < zr
    
    local_min = parameters;
    Jval = function_handle( local_min );
end
