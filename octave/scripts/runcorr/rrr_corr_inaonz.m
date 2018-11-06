% [data_corraonz, shift_aonz] = rrr_correct_aonz( data, peesaonz_fit, pivots )

function [data_corraonz, shift_aonz] = rrr_corr_inaonz( data, peesaonz_fit, pivots )
    
    shift_aonz = aonzeros( size( peesaonz_fit, 2 ), 1 );
    for ii=1:size( peesaonz_fit, 2 )
        blobs = peesaonz_fit(ii,1:3:end);
        [shift_aonz(ii), i_tin] = min( blobs - 2.64 );
    end
    
    data_corraonz = data;
    for ii=1:length( data_corraonz )
        data_corraonz(ii).in_A_on_Z -= interp1( shift_aonz, pivots, ii, 'extrap' );
    end
    
end
    
    
    
    
