%this function does the same as hybridizer.m, only ohoefully faster.
%WARNING: the interface is a little different
%
% function spc_model = hybridizer_fast( pees, hspcs, [hbkg] )
%
%parameters:
% pees : the wheighting of the spectra.
% hspc : the HISTOGRAM of the spectra, arranged in a nice matrix.
% hbkg : optionally, the HISTOGRAM of the background.
%        This will always be scaled to 1. For now.
%returns:
% spc_model : the result histogram
%
%NOTE: this function is essentially a sum, might not even be worth it

function spc_model = hybridizer_fast( pees, hspc, hbkg )
    if length( pees ) ~= size( hspc, 1 )
        error( 'pees and hspc sizes do not match up' );
    end
    
    pees = pees(:)';
    spc_model = pees*hspc;
    
    if nargin == 3
        spc_model += hbkg;
    end
end
