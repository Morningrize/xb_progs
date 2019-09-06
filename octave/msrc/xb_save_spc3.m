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
% -- ellipse: an elliptical cut, which is also optional.
% -- comments: a number of strings to put in the header

function xb_save_spc3( fname, hst, bins, herr, ell, varargin )
    if nargin == 3
        hist_errors = sqrt( hst );
    end
    
    of = fopen( [fname,'.dat'], 'w' );
    if of == -1
        error( 'Cannot open or create the file.' );
    end
    __write_header( of, varargin );
    __write_data3( of, hst, bins );
    if exist( 'ell' ); __write_cut3( fname, ell ); end
    
    fflush( of );
    fclose( of );
    
    __write_plg3( fname, hst, bins ); 
end

function __write_header( of, args )
    fprintf( of, '# File generated with XB progs, for gnuplot. YXZ row blocks.\n\n' );
    for ii=1:numel( args )
        if ischar( args{ii} )
           fprintf( of, '# %s\n', args{ii} );
        end
        fprintf( of, '\n' );
    end
end

function __write_data3( of, hst, bins )
    ybins = bins{2};
    xbins = bins{1};

    for ii=1:numel( ybins )
        fprintf( of, '%f %f %f\n', [ones( numel( xbins ), 1 )*ybins(ii), xbins(:), hst(ii,:)(:)]' );
        fprintf( of, '\n' );
    end
end

function __write_cut3( fname, ell )
    t = linspace( 0, 2*pi, 100 ); t = t(:)';
    rotM = [ cos( ell.rot ), -sin( ell.rot ); sin( ell.rot ), cos( ell.rot ) ];
    scM = [ell.a, 0; 0, ell.b];
    xy = rotM*scM*[cos(t); sin(t)] + ell.ctr(:);

    of = fopen( [fname,'-cut.dat'], 'w' );
    if of == -1; error( 'Cannot open or create the cut file.' ); end
    fprintf( of, '# File generated with XB progs, for gnuplot, XY for ellispe cut.\n\n' );
    fprintf( of, '%f %f\n', xy );
    fflush( of );
    fclose( of );
end

function __write_plg3( fname, hst, bins, titl, datalabel, xlabl, ylabl )
    %if exist( [fname,'.plg'], 'file' ); return; end

    boffx = 0.5*(bins{1}(2)-bins{1}(1));
    boffy = 0.5*(bins{2}(2)-bins{2}(1));
    if ~exist( 'titl' ); titl = fname; end
    if ~exist( 'xlabl' ); xlabl = 'X'; end
    if ~exist( 'ylabl' ); ylabl = 'Y'; end
    if ~exist( 'datalabel' ); datalabel = 'Spc'; end

    plg = fopen( [fname,'.plg'], 'w' );
    
    fprintf( plg, '# !/usr/bin/gnuplot -p -c\n\n' );
    fprintf( plg, '#This file is a starter plot, automatically generated. Edit at will!\n\n' );

    fprintf( plg, '#Some preliminaries\n' );
    fprintf( plg, 'set grid\n' );
    fprintf( plg, 'set title "%s"\n', titl );
    fprintf( plg, 'set xlabel "%s"\n', xlabl );
    fprintf( plg, 'set ylabel "%s"\n', ylabl );
    fprintf( plg, 'set xrange [%f:%f]\n', bins{1}(1)-boffx, bins{1}(end)-boffx );
    fprintf( plg, 'set yrange [%f:%f]\n', bins{2}(1)-boffy, bins{2}(end)-boffy );
    fprintf( plg, 'set cbrange [%d:%d]\n', min( min( hst) ), max( max( hst ) ) );
    fprintf( plg, 'set format cb "%%.0f"\n' );
    fprintf( plg, 'set cbtics scale 0\n' );
    fprintf( plg, 'set palette defined ( 0 "white", 17 "dark-blue", 33 "blue", 50 "yellow", 66 "orange", 100 "dark-red" )\n' );
    fprintf( plg, 'set view map\n\n' );

    fprintf( plg, '#The business end\n' );
    fprintf( plg, 'set terminal qt size 1920,1600 enhanced font "Verdana,24"\n' );
    fprintf( plg, 'plot "%s" u 2:1:3 with image title "%s", \\\n"%s" w lines lc "black" lw 2 title "cut"\n\n', [fname,'.dat'], datalabel, [fname,'-cut.dat'] );

    fprintf( plg, 'set terminal png size 1920,1600 enhanced font "Verdana,24"\n' );
    fprintf( plg, 'set output "%s"\nreplot\n\n', [fname,'.png'] );

    fprintf( plg, 'set terminal svg size 1920,1600 fname "Verdana" fsize 24\n' );
    fprintf( plg, 'set output "%s"\nreplot\n\n', [fname,'.svg'] );

    fflush( plg );
    fclose( plg );

    system( ['chmod +x ',fname,'.plg'] );
end
