%this tiny utility prints a spectrum in gnuplot-friendly format
%
% xb_save_spc( fname, hst, bins, hist_errors )
%
%parameters:
% -- fname: the name of the file to save into
% -- hst: the counts of the histogram
% -- bins: the bins of the histogram
% -- hist_errors: the error bars (it's optional, if not specified they are calculated)

function xb_save_spc( fname, hst, bins, herr )
    if nargin == 3
        hist_errors = sqrt( hst );
    end
    
    of = fopen( fname, 'w' );
    fprintf( of, '%f %f %f\n', [bins(:), hst(:), herr(:)]' );
    fflush( of );
    fclose( of );
end
    
