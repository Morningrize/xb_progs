%This utility produces a conversion factor between raw counts and reaction probability
%Basically, it handles in Tpat frequencies.
%-
% [scalers, indexes] = xb_tpat2scaler( target_dataset, [tpat_in], [tpat_out] )
%
%arguments:
% -- target_dataset: the dataset for which the spectra will be done.
% -- [tpat_in] : an array of BINARY MASKS representing the Tpats to be SELECTED
% -- [tpat_out] : an array of BINARY MASKS representing the Tpats to be DISCARDED
%returns:
% -- scalers: it's a structure with scaling information:
%             .a : the constant to scale the spectra of target_dataset( indexes.in );
%             .aa : same as above, except target_dataset( indexes.inNout );
%             .p : the product of all the "each", which is the theoretical limit
%                  if all of the downscales were independent.
%             .each : the same, for each single tpat.
%NOTE: these are effectively the downscale factons.
% -- indexes: it's a structure containing various indicizations of the dataset:
%             .in : indexes where tpat_in is exactly matched
%             .out : indexes where tpat_out is exactly excluded
%             .inNout : ins that are not out.
%             .also : indexes where tpat_in are also matched
%             .each : an array of indexes where each single tpat is matched, separately.
%
%NOTE: tpat names and flags:
%      POS_NOT_ROLU 0x0001
%      PNR_PUP      0x0003 //0x3 is not a typo, it's the actual flag
%      FRAG         0x0004
%      LAND         0x0008
%      FRAG_XB_SUMF 0x0010
%      FRAG_XB_SUM  0x0020
%      FRAG_XB_OR   0x0040
%      PIX          0x0080
%      LAND_COSM    0x0100
%      TFW_COSM     0x0200
%      NTF_COSM     0x0400
%      XB_COSM      0x0800
%      XB_SUM_OFFSP 0x1000
%      PIX_OFFSP    0x2000

function [scalers, idx] = xb_tpat2scaler( target_dataset, tpat_in, tpat_out )
    if nargin == 1
        scalers.a = 1; idx = [1:numel( target_dataset )];
        return
    end
    
    try
        tpat = [target_dataset.tpat];
    catch
        error( 'The target dataset has no Tpat in it.' );
    end
    
    tmask = 0;
    scalers.each = ones( size( tpat_in ) );
    idx.each = cell( size( tpat_in ) );
    for ii=1:numel( tpat_in )
        tmask = bitor( tpat_in(ii), tmask );
        idx.each(ii) = find( bitand( tpat, tpat_in(ii) ) );
        scalers.each(ii) = numel( target_dataset )/numel( idx.each{ii} );
    end
    scalers.p = prod( scalers.each );
    
    idx.in = find( bitand( tpat, tmask ) == tmask );
    scalers.a = numel( target_dataset )/numel( idx.in );
    idx.also = find( bitand( tpat, tmask ) );
    notmask = 0;
    
    if nargin == 3
        for ii=1:numel( tpat_out )
            notmask = bitor( tpat_out(ii), notmask );
        end
        idx.out = find( not( bitand( tpat, notmask ) == notmask ) );
        idx.inNout = intersect( idx.out, idx.in );
        scalers.aa = numel( target_dataset )/numel( idx.inNout );
    end
end
