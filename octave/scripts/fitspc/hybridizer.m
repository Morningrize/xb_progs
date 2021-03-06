%This function creates a mixture of events from some specra
%With preponderance fixed by a parameter.
%
% hybrid_spectrum = hybridizer( weights, spectra, nb_events, binZ, [background] )
%
%parameters:
% weights : the scale parameters for each spectra (must be at least as many)
% spectra : the spectra
% nb_events : the number of events the end product must have
% binZ : the binnage to actually create the spectrum.
% background : optionally, the HISTOGRAM of a background. This won't participate into the fit.
%returns:
% hybrid_spectrum : a dataset composed by the spectra, weighted,
%                   and optionally the background

function spc_model = hybridizer( pees, spectra, nbe, binZ, hbkg )
    if numel( pees ) ~= numel( spectra )
        error( 'pees and spectra must have the same numel' );
    end
    
    spc_model = [];
    for ii=1:numel( spectra )
        nn = min( round( nbe*pees(ii) ), numel( spectra{ii} ) );
        hyspc = [spc_model; ...
                 spectra{ii}(randperm( numel( spectra{ii} ), nn ))];
    end
    
    try
        nrg = xb_cluster_nrg( hyspc );
    catch
        nrg = xb_data_nrg( hyspc );
    end
    
    spc_model = hist( nrg, binZ );
    if exist( 'hbkg', 'var' ); spc_model += hbkg;
end
