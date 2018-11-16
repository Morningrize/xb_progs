% this tiny tool finds the relevant ellipse (and maybe also cuts a dataset)
%
% ell = xb_find_isoID( some_ellispes, eleA, eleZ )
% [ell, data_cut] = xb_find_isoID( some_ellispes, eleA, eleZ, data )
%
% parameters:
% -- some_ellispes: a structure array of ellispes
% -- data: OPTIONAL and switches the return types, an XB structure data set
% -- eleA: the mass number to target
% -- eleZ: the charge number to target
% returns:
% -- ell: the relevant ellispes
% -- data_cut: the elements of the cut inside the relevant ellispe (only if data is given)

function [ell, data_cut] = xb_find_isoID( some_ellispes, eleA, eleZ, data )
    if ~isstruct( some_ellispes ) || ~isfield( some_ellispes, 'ctr' ) || ...
       ~isfield( some_ellispes, 'a' ) || ~isfield( some_ellispes, 'b' )
        error( 'Not ellispses where ellispes needed.' );
    end
        
    if nargin == 3
        data = 0;
    end
    
    some_ellispes = some_ellispes;
    ctrZ = [some_ellispes.ctr];
    if size( ctrZ ) ~= [2, length( some_ellispes ) ]
        ctrZ = reshape( ctrZ, 2, length( ctrZ )/2 );
    end
    target_ctr = [eleZ, eleA/eleZ];
    
    dist = sqrt( sum( (ctrZ - target_ctr').^2 ) );
    [m, im] = min( dist );
    
    ell = some_ellispes(im);
    if nargout == 1
        return;
    end
    if ~isstruct( data ) || ~isfield( data, 'in_Z' ) ~isfield( data, 'in_A_on_Z' )
        error( 'Expected and XB data structure, got something else.' );
    end
    
    idx_cut = xb_is_in_ellipse( double( [data.in_Z] ), double( [data.in_A_on_Z] ), ell );
    data_cut = data( idx_cut );
end
