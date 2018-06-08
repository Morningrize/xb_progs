%This version of the main script has limited functionality: just
%the thresholds will be done.
%
% crysthr( run_file_list, ofile_name, [interactive] )
%
% -- run_file_list: a character matrix where each line is the name of a file
% -- ofile_name: a string used as the output file name
% -- interactive: flag to ask for user opinion

function crysthr( run_file_list, ofile_name, interactive )
	%plan:
	%0) read the file(s)
	%   NOTE: sadly the drone only deals in clusters as for now, and that gets in my way.
	%         So no drone here. sigh. This also means that the number of events that can be
	%         handled is limited.
	% 1) iteratively extract the thresholds from each crystals
	% 3) print it
	
	if ~ischar( run_file_list )
		error( 'File names must be a charachter matrix!' );
	end
	if ~ischar( ofile_name )
		error( 'Output file name must be a string!' );
	end
	
	global settings;
	settings.ax_lb = 0;
	settings.ax_ub = 3e3;
	settings.crys_nb = 0;
	settings.c_trg = 10;
	if nargin == 3
		settings.interactive = interactive;
	else
		settings.interactive = true;
	end
	
	%collect the results
	lore = struct( 'hst', cell( 1, 162 ), ...
	               'herr', cell( 1, 162 ), ...
	               'binZ', cell( 1, 162 ), ...
	               'cutoff', cell( 1, 162 ) );
	
	nf = size( run_file_list, 1 );
	data = [];
	for ii=1:nf
		data = [data, xb_load_data( deblank( run_file_list(ii,:) ) )];
	end
	
	for ii=1:162
		oh = @(p) p == ii;
		settings.crys_nb = ii;
		dcut = xb_data_cut_on_field( data, oh, 'i' );
		nrg = xb_data_nrg( dcut ); nrg = nrg( find( nrg ) );
		[lore(ii).hst, lore(ii).binZ] = hist( nrg, max( nrg )/10 );
		lore(ii).herr = sqrt( lore(ii).hst );
		
		lore(ii).cutoff = cc_do_cutoff( [lore(ii).binZ; lore(ii).hst] );
		clear dcut;
	end
	
	save( [ofile_name,'__lore'], 'lore' );
	cc_print( ofile_name, lore );
end
	
