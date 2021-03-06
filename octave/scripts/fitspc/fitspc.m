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
nrgizer = @xb_cluster_nrg;
hor_field = 'centroid_id';
binZ = [0:50:1.2e4]; %standard binnage
extremes = [1,241];
fin = {};
minopts = {};
do_fast = false;
fixed_mt = true;
fit_engine = @fitter;
do_amp = 1;
fbkg = '/home/gatto/PhD_local/xb_data/bkg-sweep-narrow/sn-132/output_files/r9xx_native_hybrid_bkg.kb.dp.xb';

function val = __check_arg( option, args, ii )
    try
        val = args{ii+1};
    catch
        error( ['option ',option,' needs an argumen.'] );
        val = 'NULL';
    end
end

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
        case { '-R', '--resume' }
            oldfile = __check_arg( '-R', args, ii );
            spc_pees = [load( oldfile, 'spc_pees' ).spc_pees];
        case { '-o', '--output-file' }
            fout = __check_arg( '-o', args, ii );
        case { '-t', '--type' }
            loader = __check_arg( '-t', args, ii );
            switch( loader )
                case { 'data', 'Data', 'xb_data' }
                    loader = @xb_load_data;
                    cutter = @xb_data_cut_on_field;
                    nrgizer = @xb_data_nrg;
                    hor_field = 'i';
                case { 'clusters', 'klz', 'cluster' }
                    %we are good already
                    ;
                case { 'adata', 'arbitrary' }
                    loader = @xb_load_adata;
                    cutter = @xb_data_cut_on_field;
                    nrgizer = @xb_data_nrg;
                    hor_field = 'i';
                otherwise
                    error( [loader,' is an unsupported format'] );
            end
        case { '-B', '--use-background' }
            fbkg = __check_arg( '-B', args, ii );
        case { '-nB', '--no-background' }
            clear fbkg;
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
        case { '-F', '--fast' }
            do_fast = true;
	case { '-f', '--fastish' }
            do_fastish = true;
        case { '-m.lr', '--minimizer.learning-rate' }
            minopts = [minopts, 'lr', str2num( __check_arg( '-m.lr', args, ii ) )];
        case { '-m.zr', '--minimizer.zero-is' }
            minopts = [minopts, 'z', str2num( __check_arg( '-m.z', args, ii ) )];
        case { '-m.M', '--minimizer.max-iter' }
            minopts = [minopts, 'M', str2num( __check_arg( '-m.M', args, ii ) )];
        case { '-M', '--manual' }
            fit_engine = @manual_fitter;
        case { '-nMT', '--no-fixed-MT' }
            fixed_mt = false;
        case { '-A', '--data-amp' }
            do_amp = str2num( __check_arg( '-a', args, ii ) );
    end
end

if ~fixed_mt && ~do_fast && ~do_fastish
    error( 'MT must be fixed to resample the event spectra.' )
end

%Assuming the data file is the last of the list.
data = loader( fin{end} );
disp( 'Data on board!' );
fin = fin(1:end-1);
[h_data, ~, herr_data] = xb_make_spc_ffb( data, binZ );
if do_amp ~= 1
    for ii=1:3
        [h_data(ii), herr_data(ii)] = xb_spcamp( h_data{ii}, binZ, do_amp );
    end
    disp( ['Data amplified by ',num2str(do_amp)] );
end

spectra = {};
for ii=1:numel( fin )
    spectra(ii) = loader( fin{ii} );
    %let's try to save a BIT of RAM
    if numel( spectra{ii} ) > numel( data )*do_amp
        spectra{ii} = spectra{ii}( randperm( numel( spectra{ii} ), numel( data )*do_amp ) );
    end
    disp( ['Spectrum in file "',fin{ii},'" loaded and cut to measure'] );
end

if exist( 'fbkg', 'var' )
    bkg = loader( fbkg );
    [hbkg, ~, herrbkg] = xb_make_spc_ffb( bkg, binZ );
    if numel( bkg ) < numel( data )*do_amp
        for ii=1:3
            [hbkg(ii), herrbkg(ii)] = xb_spcamp( hbkg{ii}, binZ, numel(data)/numel(bkg)*do_amp );
        end
    end
    disp( 'Background loaded.' );
end

%build the empty target model
%NOTE: this is valid only for the front of the CB.
%      full and back may become available later.
if fixed_mt
    load /home/gatto/PhD_local/xb_data/xb/empty/9xx_AR/v01/empty-target-functional
    mtf = mt_scaler_2s133s( numel( data )*do_amp )*mt_model( binZ, pees_front );
    disp( 'MT target model loaded.' );
else
    spectra(end+1) = loader( '/home/gatto/PhD_local/xb_data/xb/empty/9xx_AR/v01/r9xx_mt_132in_2s133s.kb.dp.xb' );
    disp( 'MT added to the spectra.' );
end

icbf = [xb_ball_neigh( 81, 5 ).i];
ohf = @(p) xb_op_cbi( p, icbf );
for ii=1:numel( spectra )
    spectra(ii) = cutter( spectra{ii}, ohf, hor_field );
end
disp( 'Spectra cut to front XB.' );

%make the bkg
if exist( 'hbkg', 'var' ) && fixed_mt
    hbkg_f = hbkg{2} + round( mtf );
elseif ~exist( 'hbkg', 'var' ) && fixed_mt
    hbkg_f = round( mtf );
elseif exist( 'hbkg', 'var' ) && ~fixed_mt
    hbkg_f = hbkg{2};
end

%TODO list:
% 0) cut all the spectra to only contain the _front_
% 0.5) if you have the background, frontize it as well.
% 1) hybridize the thing with _variable_ parameters
% 2) do the spectrum
% 3) fit 1-2 to the data (also, a part of the data if needed).
% 4) Do all this efficiently.

if ~exist( 'spc_pees', 'var' ); spc_pees = ones( 1, numel( spectra ) ); end 
if do_fast || do_fastish
    hspc = [];
    for ii=1:numel( spectra )
        nrg = nrgizer( spectra{ii} );
        hspc = [hspc; hist( nrg, binZ )];
    end
    if ~fixed_mt
        pkg load signal
        hspc(end,:) = sgolayfilt( [hspc(end,1:end-1),0], 2, 21 );
    end
    if do_fast; spc_model = @( pees ) hybridizer_fast( pees, hspc, hbkg_f );
    elseif do_fastish;
        spc_model = @( pees ) hybridizer_fastish( pees, hspc, hbkg_f, binZ );
    end
else
    spc_model = @( pees ) hybridizer( pees, spectra, numel( data ), binZ, hbkg_f );
end

if isempty( minopts )
    [spc_pees, spc_errs, chisq] = fit_engine( spc_pees, spc_model, h_data{2}, extremes, binZ );
else
    [spc_pees, spc_errs, chisq] = fit_engine( spc_pees, spc_model, h_data{2}, extremes, ...
                                               binZ, minopts );
end

%make also the graph (the extra mile bit)
spc_modelerr = @( p ) sqrt( hybridizer_fast( p.^2, hspc.^2 ) + hbkg_f );
spc_err_zeroed = zeros( size( spc_errs ) );
spc_err_zeroed( find( spc_pees ) ) = spc_errs( find( spc_pees ) );
xb_save_spc( fout, { h_data{2}, spc_model( spc_pees ) }, ...
             { binZ, binZ }, { sqrt( h_data{2} ), spc_modelerr( spc_err_zeroed ) } );

%writeout time
save( fout, 'fin', 'spc_pees', 'spc_errs', 'spc_model', 'spc_modelerr', 'chisq' );
