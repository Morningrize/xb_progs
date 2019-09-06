%This utility takes in what would be a scatter plot and does a nice 3d histogram
%It requires the statistics package in Octave. It checks for it.
%
% [hst, binz, herr] = xb_make_spc3( XYdata, bin, div )
%
%parameters:
% -- XYdata: a many-by2 matrix, first column X, second column Y.
% -- bin: the description of the binnage: a rectangle, or the same as hist3.
% -- div: the number of divisions if bin is a rectangle.
%returns:
% -- hst: the big fat matrix with all the counts in it
% -- binz: the centres of the bins, as returned by hist3
% -- herr: the big fat matrix of the errors (useless for now, maybe in the future)

function [hst, binz, herr] = xb_make_spc3( pairs, bin, div )
    if nargin == 2
        div = 512;
    end

    try
        pkg load statistics
    catch
        warning( 'Attempting to install the statistics package' );
        pkg install -forge statistics
        pkg load statistics
    end

    if size( pairs, 2 ) ~= 2 && size( pairs, 1 ) == 2
        pairs = pairs';
    else
        error( 'pairs must be a 2-by-many or many-by-2 matrix.' );
    end
    
    if nargout == 0
        if isscalar( bin )
            hist3( pairs, [bin, bin] );
        elseif isvector( bin ) && numel( bin ) == 4
            hist3( pairs, { linspace( bin(1), bin(2), div ), linspace( bin(3), bin(4), div ) } );
        elseif iscell( bin )
            hist3( pairs, bin );
        end
    else
        if isscalar( bin )
            [hst, binz] = hist3( pairs, [bin, bin] );
        elseif isvector( bin ) && numel( bin ) == 4
            [hst, binz] = hist3( pairs, { linspace( bin(1), bin(2), div ), linspace( bin(3), bin(4), div ) } );
        elseif iscell( bin )
            [hst, binz] = hist3( pairs, bin );
        end
        hst = hst'; %this appears to be necessary to make things the right way up
        herr = sqrt( hst );
    end
end
