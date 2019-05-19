%this function does the same as hybridizer.m, only ohoefully faster.
%WARNING: the interface is a little different
%
% function spc_model = hybridizer_fast( pees, hspcs, [hbkg] )
%
%parameters:
% pees : the wheighting of the spectra.
% hspc : the HISTOGRAM of the spectra, in a cell array.
% hbkg : optionally, the HISTOGRAM of the background.
%        This will always be scaled to 1. For now.
%returns:
% spc_model : the result histogram
%
%NOTE: this function is essentially a sum, might not even be worth it

function spc_model = hybridizer_fast( pees, hspc, hbkg )
    if numel( pees ) ~= numel( hspc )
        error( 'pees and hspc must have the same numel' );
    end
    
    spc_model = zeros( size( hspc{1} ) );
    for ii=1:numel( hspc )
        spc_model += round( pees(ii)*hspc{ii} );
    end
    
    if nargin == 3
        spc_model += hbkg;
    end
end
