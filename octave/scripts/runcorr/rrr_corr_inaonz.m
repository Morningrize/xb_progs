%This function corrects in_A_on_Z to center the whole thing on the "target" position
%If  "target" is omitted, it defaults to 2.64
%
% [data_corraonz, shift_aonz] = rrr_corr_inaonz( data, peesaonz_fit, pivots, target )

function [data, shift_aonz] = rrr_corr_inaonz( data, peesaonz_fit, pivots, target )
    if nargin == 3
        target = 2.64;
    end
    
    shift_aonz = zeros( size( peesaonz_fit, 2 ), 1 );
    for ii=1:size( peesaonz_fit, 2 )
        blobs = peesaonz_fit{ii};
        blobs = blobs(2:3:end);
        [shift_aonz(ii), i_tin] = min( abs( blobs - target ) );
    end
    
    coeff = xb_OLS( pivots, shift_aonz, 2 );
    shifun = @( p ) coeff(1) + coeff(2)*p + coeff(3)*p.^2;
    
    inaonz_shift = shifun( [1:length(data)] );
    for ii=1:length( data )
        data(ii).in_A_on_Z += inaonz_shift(ii);
    end
    
end
    
    
    
    
