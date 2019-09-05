%this tiny utility prints a spectrum in gnuplot-friendly format, the 3D edition!!!
%And will also generate a .plg file to build the spectrum (you should then customize it in there)
%
% xb_save_spc3( fname, hst, bins, hist_errors, comments )
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
    
    of = fopen( [fname,'.dat'], 'w' );
    if of == -1
        error( 'Cannot open or create the file.' );
    end
    __write_header( of, varargin );

    dlmwrite( of, hst, '    ', 0, 0 );
    fflush( of );
    fclose( of );
    
    __write_plg3( fname, [max( max(hst) ), bins{1}(1), bins{1}(end), bins{2}(1), bins{2}(end)] ); 
end

function __write_header( of, args )
    fprintf( of, '# File generated with XB progs, for gnuplot. XY heatmap of Z.\n' );
    for ii=1:numel( args )
        if ischar( args{ii} )
           fprintf( of, '# %s\n', args{ii} );
        end
        fprintf( of, '\n' );
    end
end

function __write_plg3( fname, dspec, titl, datalabel, xlabl, ylabl )
    if exist( [fname,'.plg'], 'file' ); return; end

    if ~exist( 'titl' ); titl = fname; end
    if ~exist( 'xlabl' ); xlabl = 'X'; end
    if ~exist( 'ylabl' ); ylabl = 'Y'; end
    if ~exist( 'datalabel' ); datalabel = 'Spc'; end

    plg = fopen( [fname,'.plg'], 'w' );
    
    fprintf( plg, '#! /usr/bin/gnuplot -p -c\n\n' );
    fprintf( plg, '#This file is a starter plot, automatically generated. Edit at will!\n\n' );

    fprintf( plg, '#Some preliminaries\n' );
    fprintf( plg, 'set grid\n' );
    fprintf( plg, 'set title "%s"\n', titl );
    fprintf( plg, 'set xlabel "%s"\n', xlabl );
    fprintf( plg, 'set ylabel "%s"\n', ylabl );
    fprintf( plg, 'set xrange [%f:%f]\n', dspec(2), dpsec(3) );
    fprintf( plg, 'set yrange [%f:%f]\n', dspec(4), dspec(5) );
    fprintf( plg, 'set cbrange [0:%d]\n', dspec(1) );
    fprintf( plg, 'set cb format "%%d"\n' );
    fprintf( plg, 'set cbtics' );
    fprintf( plg, 'set palette defined ( 0 "white", 17 "dark blue", 33 "blue", 50 "yellow", 66 "red", 100 "dark red" )\n' );
    fprintf( plg, 'set view map\n' );

    fprintf( plg, '#The business end\n' );
    fprintf( plg, 'set terminal qt size 1920,1600 enhanced font "Verdana,24"\n' );
    fprintf( plg, 'splot ' );
    fprintf( plg, '"%s" matrix title "%s" w image\n\n', [fname,'.dat'], datalabel );

    fprintf( plg, 'set terminal png size 1920,1600 enhanced font "Verdana,24"\n' );
    fprintf( plg, 'set output "%s"\nreplot\n\n', [fname,'.png'] );

    fprintf( plg, 'set terminal svg size 1920,1600 fname "Verdana" fsize 24\n' );
    fprintf( plg, 'set output "%s"\nreplot\n\n', [fname,'.svg'] );

    fflush( plg );
    fclose( plg );

    system( ['chmod +x ',fname,'.plg'] );
end
