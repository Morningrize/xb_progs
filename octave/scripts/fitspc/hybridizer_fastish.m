%This tool is yet another hybridizer, hopefully faster than the resampler and
%not too muc slower than the fast one --well, it will behave like the fast one
%in some sense
%
% spc_model = hybridiszer_fastish( pees, hspcs, [hbkg, binZ] )
%
%parameters:
% pees : the wheighting of the spectra.
% hspc : the HISTOGRAM of the spectra, arranged in a nice matrix.
% hbkg : optionally, the HISTOGRAM of the background.
% binZ : optionally, the COMMON binnage --defaults at [0:50:12e3]
%returns:
% spc_model : the result histogram

function spc_model = hybridiszer_fastish( pees, hspc, hbkg, binZ )
    if length( pees ) ~= size( hspc, 1 )
        error( 'pees and hspc size do not match up' );
    end

    if ~exist( 'binZ', 'var' ); binZ = [0:50:12e3]; end

    for ii=1:numel( pees )
        hspc(ii,:) = xb_spcamp( hspc(ii,:), binZ, pees(ii) );
    end

    spc_model = sum( hspc, 1 );
    if exist( 'hbkg', 'var' ); spc_model += hbkg; end
end
