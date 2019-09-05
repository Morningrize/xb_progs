%this tiny utility prints a spectrum in gnuplot-friendly format
%And will also generate a .plg file to build the spectrum (you should then customize it in there)
%
% xb_save_spc( fname, hst, bins, hist_errors, comments )
%
%parameters:
% -- fname: the name of the file to save into. Possibly without extension!
% -- hst: the counts of the histogram
% -- bins: the bins of the histogram
% -- hist_errors: the error bars (it's optional, if not specified they are calculated)
% -- comments: a number of strings to put in the header

function xb_save_spc( fname, hst, bins, herr, varargin )
    if nargin == 3
        hist_errors = sqrt( hst );
    end
    
    if iscell( bins )
        fnames = cell( size( bins ) );
        nb_dat = numel( bins );
        for ii=1:nb_dat
            fnames(ii) = [fname,'-',num2str(ii)];
            __write_data( fnames{ii}, bins{ii}, hst{ii}, herr{ii} );
        end
    else
        nb_dat = 1;
        __write_data( [fname,'-1'], bins, hst, herr );
    end

    __write_plg( fname, nb_dat ); 
end

function __write_header( of, args )
    fprintf( of, '# File generated with XB progs, for gnuplot. | Bin | Count | sqrt( Count ) |\n' );
    for ii=1:numel( args )
        if ischar( args{ii} )
           fprintf( of, '# %s\n', args{ii} );
        end
        fprintf( of, '\n' );
    end
end

function __write_data( fname, bins, hst, herr, varargin )
    of = fopen( [fname,'.dat'], 'w' );
    if of == -1
        error( 'Cannot open or create the file.' );
    end
    __write_header( of, varargin );

    fprintf( of, '%f %f %f\n', [bins(:), hst(:), herr(:)]' );
    fflush( of );
    fclose( of );
end

function __write_plg( fname, nb_files, titl, datalabel, xlabl, ylabl )
    if exist( [fname,'.plg'], 'file' ); return; end

    if ~exist( 'titl' ); titl = fname; end
    if ~exist( 'xlabl' ); xlabl = 'KeV'; end
    if ~exist( 'ylabl' ); ylabl = '#/Kev'; end
    if ~exist( 'datalabel' ); datalabel = 'Spc'; end

    plg = fopen( [fname,'.plg'], 'w' );
    
    fprintf( plg, '#! /usr/bin/gnuplot -p -c\n\n' );
    fprintf( plg, '#This file is a starter plot, automatically generated. Edit at will!\n\n' );

    fprintf( plg, '#Some preliminaries\n' );
    fprintf( plg, 'set grid\n' );
    fprintf( plg, 'set title "%s"\n', titl );
    fprintf( plg, 'set xlabel "%s"\n', xlabl );
    fprintf( plg, 'set ylabel "%s"\n', ylabl );
    fprintf( plg, 'set logscale y\n' );
    fprintf( plg, 'set bars fullwidth\n\n' );

    fprintf( plg, '#The business end\n' );
    fprintf( plg, 'set terminal qt size 1920,1600 enhanced font "Verdana,24"\n' );
    fprintf( plg, 'plot ' );
    for ii=1:nb_files
        fprintf( plg, '"%s" u 1:2 title "%s" w histeps lw 2 lc %d, "" u 1:2:3 notitle w errorbars ps 0 lc %d', [fname,'-',num2str(ii),'.dat'], datalabel, ii, ii );
        if ii < nb_files; fprintf( plg, ', \\\n' ); end
    end
    fprintf( plg, '\n\n' );

    fprintf( plg, 'set terminal png size 1920,1600 enhanced font "Verdana,24"\n' );
    fprintf( plg, 'set output "%s"\nreplot\n\n', [fname,'.png'] );

    fprintf( plg, 'set terminal svg size 1920,1600 fname "Verdana" fsize 24\n' );
    fprintf( plg, 'set output "%s"\nreplot\n\n', [fname,'.svg'] );

    fflush( plg );
    fclose( plg );

    system( ['chmod +x ',fname,'.plg'] );
end
