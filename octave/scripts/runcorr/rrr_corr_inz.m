% This function corrects a dataset interpolating on where the isotope blob for tin (z=50)
% is supposed to be
%
% [data_corrz, shift_z] = rrr_correct_z( data, peesZ_fit, pivots )

function [data_corrz, shift_z] = rrr_corr_inz( data, peesZ_fit, pivots )
    
    shift_z = zeros( length( peesZ_fit ), 1 );
    for ii=1:length( peesZ_fit );
        blobs = peesZ_fit{ii};
        blobs = blobs(2:3:end);
        [shift_z(ii), i_tin] = min( blobs - 50 );
    end
    
     data_corrz = data;
     for ii=1:length( data_corrz )
         data_corrz(ii).in_Z -= interp1( shift_z, pivots, ii, 'extrap' );
     end
    
end
    
    
