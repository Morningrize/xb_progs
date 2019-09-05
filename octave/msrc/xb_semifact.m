%this utility calculates (hopefully, not too slowly) the semifactorial
%of the argument.
%
%NOTE: n!! is freaking not (n!)!, but this function. Just to keep confusione high.
%
% val = xb_semifact( x )
%
%parameters:
% x : the number to semifactorialize (it's ought to be an integer)
%returns:
% val : the result of the operation.

function val = xb_semifact( x )
    [i_even] = find( mod( x, 2 ) == 0 );
    [i_odd] = find( mod( x, 2 ) );
    
    val = zeros( size( x ) );
    k = x(i_even)/2;
    val(i_even) = 2.^k.*factorial( k );
    k = (x(i_odd)+1)/2;
    val(i_odd) = factorial( 2*k-1 )./factorial( k-1 ).*2.^(1-k);
end
    
