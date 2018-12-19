%This main gets a background simulation and produces a file that cointains
%the average illumination for every crystal.
%
% crysnoise( noise_file_list, ofile_name, options )
%
%Parameters:
% -- noise_file_list: it's a list of one or more noise files
% -- ofile_name: it's the name of the output file
% -- options: many many exciting options....
%Returns: nothing, writes the file

function crysnoise( noise_file_list, ofile_name, varargin )
    %plan:
    %0) read the noise files --I'm assuming these are data files, before clustering.
    %                          The reason being: atomic bkg isn't very energetic and should
    %                          be contained in one single crystal most of the time.
    %1) extract the energies for every single crystal
    %2) Do an average and shove it into a file
    
    if ~ischar( noise_file_list )
		error( 'File names must be a charachter matrix!' );
	end
	if ~ischar( ofile_name )
		error( 'Output file name must be a string!' );
	end
	
	%collect the results
	lore = struct( 'hst', cell( 1, 162 ), ...
	               'herr', cell( 1, 162 ), ...
	               'binZ', cell( 1, 162 ), ...
	               'avg_illum', cell( 1, 162 ), ...
                   'avg_illum_phit', cell( 1, 162 )
                 );
	
	nf = size( noise_file_list, 1 );
	data = [];
	for ii=1:nf
		data = [data, xb_load_data( deblank( noise_file_list(ii,:) ) )];
	end
	
	for ii=1:162
		oh = @(p) p == ii;
		settings.crys_nb = ii;
		dcut = xb_data_cut_on_field( data, oh, 'i' );
		nrg = xb_data_nrg( dcut );
        try
            [lore(ii).hst, lore(ii).binZ] = hist( nrg, max( nrg )/50 );
            lore(ii).herr = sqrt( lore(ii).hst );
            %NOTE: pdist = h/sum( h ); phit = sum( h )/numel( data );
            %      ---> pdist*phit = h/numel( data );
            lore(ii).avg_illum_phit = sum( (lore(ii).hst.*lore(ii).binZ) )/numel( data );
            lore(ii).avg_illum = sum( (lore(ii).hst.*lore(ii).binZ) )/sum( lore(ii).hst );
        catch
            lore(ii).hst = []; lore(ii).binz = []; lore(ii).herr = [];
            lore(ii).avg_illum = 0;
        end
		
		clear dcut;
	end
	
	save( [ofile_name,'__lore'], 'lore' );
    fields = { 'avg_illum', 'avg_illum_phit' };
	cc_print( ofile_name, 'w', lore, fields );
end
