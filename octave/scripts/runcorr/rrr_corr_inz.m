% This function corrects a dataset interpolating on where the isotope blob for tin (z=50)
% is supposed to be
%
% [data_corrz, shift_z] = rrr_correct_z( data, peesZ_fit, pivots )

function [data, shift_z] = rrr_corr_inz( data, peesZ_fit, pivots )
    
    shift_z = zeros( length( peesZ_fit ), 1 );
    for ii=1:length( peesZ_fit );
        blobs = peesZ_fit{ii};
        blobs = blobs(2:3:end);
        [shift_z(ii), i_tin] = min( abs( blobs - 50 ) );
    end
    
    coeff = xb_OLS( pivots, shift_z, 2 );
    shifun = @( p ) coeff(1) + coeff(2)*p + coeff(3)*p.^2;

    inz_shift = shifun( [1:length(data)] );
    inz = double( [data.in_Z] ) - inz_shift;
    for ii=1:length( data )
        data(ii).in_Z = inz(ii);
    end
end
    
    
