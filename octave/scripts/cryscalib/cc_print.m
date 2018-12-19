%this function prints the calibration data into a hooman readable file
%
% file = cc_print( file{_name}, writing_mode, lore )
% -- lore: it's the structure where the results are stored (and that's also saved)
%          here an easier file is produced. Per crystal.
%returns nothing.
%outputs a file.

function file = cc_print( file, writing_mode, lore, fields )
	if ischar( file )
		file = fopen( file, writing_mode );
	end
	
	if nargin == 3
        fields = { 'cutoff', 'cal_p', 'cal_e', 'dE_E' };
    end
	
	for cc=1:162
		fprintf( file, 'Crystal number %d\n', cc );
		
        for ff=1:numel( fields )
            if isfield( lore, fields{ff} )
                fprintf( file, fields{ff} );
                fprintf( file, ' %f\n', lore(cc).( fields{ff} ) );
            end
        end
		
		fprintf( file, '\n@@@\n\n' );
	end

	fflush( file );
end
