#! /usr/bin/octave --no-gui

#Point this script at a number of files containing spectra, optionally specify the binnage
#And enjoy a completely pointless fitting --yup, that's science.
#
# ./fitspc.m [files] [result-report]

args = [ [pwd,'/fitspc.m']; argv() ];
fin = {};
fout = 'fitted_spectra.octave';
loader = @xb_load_clusterZ;
cutter = @xb_cluster_cut_on_field;
hor_field = 'centroid_id';
binZ = [0:50:1.2e4]; %standard binnage
extremes = [0,1.2e4];
fin = {};

for ii=2:numel( args )
    curr = args{ii};
    prev = args{ii-1};
    if curr(1) ~= '-' && prev(1) ~= '-'
        fin = [fin; curr];
        continue;
    end
    
    switch( curr )
        case { '-v', '--verbose' } %has no effect ATM
            verbose = true;
        case { '-r', '--reaction-probability' }
            %NOTE: this is the reaction probability for the 2+
            %      other reaction probabilities may become available
            %      in the future.
            %NOTE: this doesn't work yet.
            rp = __check_arg( '-r', args, ii );
            error( '-r needs an argumen.' );
        case { '-o', '--output-file' }
            fout = __check_arg( '-o', args, ii );
        case { '-t', '--type' }
            loader = __check_arg( '-d', args, ii );
            switch( loader )
                case { 'data', 'Data', 'xb_data' }
                    loader = @xb_load_data;
                    cutter = @xb_data_cut_on_field;
                    hor_field = 'i';
                case { 'clusters', 'klz', 'cluster' }
                    loader = @xb_load_clusterZ;
                    cutter = @xb_cluster_cut_on_field;
                    hor_field = 'centroid_id';
                case { 'adata', 'arbitrary' }
                    loader = @xb_load_adata;
                    cutter = @xb_data_cut_on_field;
                    hor_field = 'i';
                otherwise
                    error( [loader,' is an unsupported format'] );
            end
        case { '-B', '--atomic-background' }
            fbkg = __check_arg( '-B', args, ii );
        case { '-b', '--binnage' }
            try
                binZ = sscanf( args{ii+1}, '[%f:%f:%f]' )
            catch
                error( 'option -b requires an [ar:gu:ment]' );
            end
            binZ = [binZ(1):binZ(2):binZ(3)];
        case { '-x', '--extremes' }
            extremes = __check_arg( '-x', args, ii );
            extremes = sscanf( extremes, '[%f:%f]' );
        otherwise
            error( ['option ',curr,' does not exist.'] );
    end
end

%Assuming the data file is the last of the list.
data = loader( fin{end} );
if exist( 'fbkg', 'var' )
    bkg = loader( fbkg );
    bkg = bkg( randperm( numel( bkg ), numel( data ) ) );
end

spectra = {};
for ii=1:numel( fin )-1
    spectra(ii) = loader( fin{ii} );
    %let's try to save a BIT of RAM
    if numel( spectra{ii} ) > numel( data )
        spectra{ii} = spectra{ii}( randperm( numel( spectra{ii} ), numel( data ) ) );
    end
end

%build the empty target model
%NOTE: this is valid only for the front of the CB.
%      full and back may become available later.
load /home/gatto/PhD_local/xb_data/xb/empty/9xx_AR/v00/empty-target-model
mtf = gscaler( numel( data ) )*gmodel_mtf( pees_mtf, { [], binZ } );

[h_data, ~, herr_data] = xb_make_spc_ffb( data, binZ );

%TODO list:
% 0) cut all the spectra to only contain the _front_
% 0.5) if you have the background, frontize it as well.
% 1) hybridize the thing with _variable_ parameters
% 2) do the spectrum
% 3) fit 1-2 to the data (also, a part of the data if needed).
% 4) Do all this efficiently.

icbf = [xb_ball_neigh( 81, 5 ).i];
ohf = @(p) xb_op_cbi( p, icbf );

for ii=1:numel( spectra )
    spectra(ii) = cutter( spectra{ii}, ohf, hor_field );
end

spc_pees = 0.001*ones( 1, numel( spectra ) ); 
if ~exist( 'bkg', 'var' )
    spc_model = @( pees ) hybridizer( pees, spectra, numel( data ), binZ );
else
    spc_model = @( pees ) hybridizer( pees, spectra, numel( data ), binZ, bkg );
end

spc_pees = fitter( spc_pees, spc_model, h_data{2}, binZ, extremes, mtf );

%writeout time
save( fout, 'spc_pees', 'spc_model' );

%----------------------------------------------------------------------------------

function val = __check_arg( option, args, ii )
    try
        val = args{ii+1};
    catch
        error( ['option ',option,' needs an argumen.'] );
        val = 'NULL';
    end
end
