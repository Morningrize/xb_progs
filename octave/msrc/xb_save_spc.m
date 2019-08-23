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
    
    of = fopen( [fname,'.dat'], 'w' );
    if of == -1
        error( 'Cannot open or create the file.' );
    end
    args = varargin;
    if ~isempty( args )
        __write_header( of, varargin );
    end

    fprintf( of, '%f %f %f\n', [bins(:), hst(:), herr(:)]' );
    fflush( of );
    fclose( of );

    __write_plg( fname );
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

function __write_plg( fname, title, datalabel, xlabel, ylabel )
    if ~exist( 'title' ); title = fname; done
    if ~exist( 'xlabel' ); xlabel = 'KeV';
    if ~exist( 'ylabel' ); ylabel = '#/Kev';
    if ~exist( 'datalabel' ); datalabel = 'Spc';

    plg = fopen( [fname,'.plg'], 'w' );
    
    fprintf( plg, '#! /usr/bin/gnuplot -p -c\n\n' );
    fprintf( plg, '#This file is a starter plot, automatically generated. Edit at will!\n\n' );

    fprintf( plg, '#Some preliminaries\n' );
    fprintf( plg 'set grid\n' );
    fprintf( plg, 'set title "%s"\n', title );
    fprintf( plg, 'set xlabel "%s"\n', xlabel );
    fprintf( plg, 'set ylabel "%s"\n\n', ylabel );

    fprintf( plg, '#The business end\n' );
    fprintf( plg, 'set terminal qt size 1920,1600 enhanced font "Verdana,24"\n' );
    fprintf( plg, 'plot "%s" u 1:2 title "%s" w histeps lw 2, "" u 1:2:3 notitle w errorbars ps 0\n\n', [fname,'.dat'], datalabel );

    fptintf( plg, 'set terminal png size 1920,1600 enhanced font "Verdana,24"\n' );
    fprintf( plg, 'set output "%s"\nreplot\n\n\', [fname,'.png'] );

    fprintf( plg, 'set terminal svg size 1920,1600 fname "Verdana" fsize 24\n' );
    fprintf( plg, 'set output "%s"\nreplot\n\n\', [fname,'.svg'] );

    fflush( plg );
    fclose( plg );

    system( ['chmod +x ',fname,'.plg'] );
end
