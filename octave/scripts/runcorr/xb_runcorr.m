%this fat tool will correct some run files for the shift on the incoming Z
%and then return the result, ready for cuts.
%
% [data_corrected, shift_Z, shift_AonZ, shift_Z fit, shift_AonZ_fit] = xb_runcorr( data, isotope_Z, isotope_A, event_step )
%
% parameters:
% -- data: a data/track structure array (as it comes out of the box) containing one or more runs.       
% -- isotope_Z: the target charge to fit.
% -- isotope_A: the mass we're looking for.
% -- event_step: OPTIONAL, check the shift every event_step events.
% returns:
% --data_corrected: the same dataset, shifted.
% --shift_Z: a 2xmany matrix containing the shift in the first row and the error thereof in the second
% --shift_AonZ: same as above, but with the mass to charge ratio.
% --shift_{Z,AonZ}_fit: the fit parameters and their errors.
%
%The philosophy of this corrector is to have at least 30-ish point on a dataset and interpolate a correction for each single event --assuming the events are in order and thus assuming that the event ID is in ascending order.

function [data_corrected, shift_Z, shift_AonZ, ...
          peesZ_fit, errZ_fit, peesAonZ_fit, errAonZ_fit] = xb_runcorr( data, nb_pivots )
    if nargin == 1
        nb_pivots = 30;
    end

    inz = double( [data.in_Z] );
    inaonz = double( [data.in_A_on_Z] );
    
    d_length = length( data );
    corrector_halfwidth = round( d_length/(1.7*nb_pivots) ); %I want overlap
    pivots = [1:round( d_length/nb_pivots ):d_length];
    
    peesZ_fit = {};
    errZ_fit = {};

    for pp=pivots
        idx_min = max( pp - corrector_halfwidth, 1 );
        idx_max = min( pp + corrector_halfwidth, d_length );
        
        [hinz, binz] = hist( inz(idx_min:idx_max), [min(inz):0.05:max(inz)] );
        [inz_m, inz_im] = xb_multigaus_find( hinz, 'triglevel', max( hinz )/10 );
        pees_init = zeros( 3*numel( inz_m ), 1 );
        pees_init(1:3:end) = inz_m;
        pees_init(2:3:end) = binz( inz_im );
        pees_init(3:3:end) = 0.2;
        
        [pf, ~, pe] = xb_multigaus_fit( [binz; hinz], pees_init );
        peesZ_fit(end+1) =  pf;
        errZ_fit(end+1) = pe;
    end

    [data, shift_Z] = rrr_corr_inz( data, peesZ_fit, pivots );
    target_inaonz = 132/50;
    inz = double( [data.in_Z] );

    peesAonZ_fit = {};
    errAonZ_fit = {};
    
    for pp=pivots
        idx_min = max( pp - corrector_halfwidth, 1 );
        idx_max = min( pp + corrector_halfwidth, d_length );
        selection = find( inz(idx_min:idx_max) <= 50 );
    
        [hinaonz, binaonz] = hist( inaonz(find( inz < zcut )), ...
                                   [min(inaonz):0.001:max(inaonz)] );
        [inaonz_m, inaonz_im] = xb_multigaus_find( hinaonz, 'triglevel', ...
                                                max( hinaonz )/10 );
        pees_init = zeros( 3*numel( inaonz_m ), 1 );
        pees_init(1:3:end) = inaonz_m;
        pees_init(2:3:end) = binaonz( inaonz_im );
        pees_init(3:3:end) = 0.005;

        [pf, ~, pe] = xb_multigaus_fit( [binaonz; hinaonz], pees_init );
        peesAonZ_fit(end+1) =  pf;
        errAonZ_fit(end+1) = pe;
    end
    
    [data_corrected, shift_AonZ] = rrr_corr_inaonz( data, peesAonZ_fit, ...
                                                    pivots, target_inaonz );
end
