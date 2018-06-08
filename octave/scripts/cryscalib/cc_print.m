%this function prints the calibration data into a hooman readable file
%
% file = cc_print( file{_name}, writin_mode, lore )
% -- lore: it's the structure where the results are stored (and that's also saved)
%          here an easier file is produced. Per crystal.
%returns nothing.
%outputs a file.

function file = cc_print( file, writing_mode, lore )
	if ischar( file )
		file = fopen( file, writing_mode );
	end
	
	for cc=1:162
		fprintf( file, 'Crystal number %d\n', cc );
		
		if isfield( lore, 'cutoff' )
			fprintf( file, 'Cutoff' );
			fprintf( file, ' %f', lore(cc).cutoff );
			fprintf( file, '\n' );
		end
		
		if isfield( lore, 'cal_p' )
			fprintf( file, 'Calibration' );
			fprintf( file, ' %f', lore(cc).cal_p );
			fprintf( file, '\n' );
		end
		
		if isfield( lore, 'cal_e' )
			fprintf( file, 'Errors' );
			fprintf( file, ' %f', lore(cc).cal_e );
			fprintf( file, '\n' );
		end
		
		if isfield( lore, 'dE_E' )
			fprintf( file, 'dE_E' )
			fprintf( file, ' %f', lore(cc).dE_E );
			fprintf( file, '\n' );
		end
		
		fprintf( file, '\n@@@\n\n' );
	end

	fflush( file );
end
