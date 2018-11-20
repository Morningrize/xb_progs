%This tool finds the incoming isotope ID out of an XB data structure.
%It should work with all structures, from adata to track_info.
%NOTE: it's not a requirement, but it will work _much_ better if the dataset
%      has been treated with xb_run_corr.
%NOTE: this tool has limited capabilities, i.e. it must be possible to distinguish
%      isotopes and isotones in the projection on the axis. Better thing should be coming.
%
% [ blob_ellispes, peesZ, errZ, peesAonZ, errAonZ ] = xb_isoID( data, nsigma, zcut );
% [...] = xb_isoID( [inz,inaonz], ... );
%
% parameters:
% -- data: it's an XB structure array or a many-by-two matrix
% -- nsigma: optionally, do the ellispese at N sigma instead of 2 (default)
% -- zcut: optionally, 
% returns:
% -- blob_ellispes: it's a array of all the ellipses representing the blobs
% -- pees{Z,AonZ}: the fitted parameters
% -- err{Z,AonZ}: the errors thereupon

function [blob_ellispes, peesZ, errZ, peesAonZ, errAonZ] = xb_isoID( data, nsigma )
    if isstruct( data ) && ~isfield( data, 'in_Z' ) || ~isfield( data, 'in_A_on_Z' )
        error( 'Not an XB data structure or a matrix.' );
    elseif ~ismatrix( data )
        error( 'Not a matrix nor an XB structure' )
    end
    
    if isstruct( data )
        inz = double( [data.in_Z] );
        inaonz = double( [data.in_A_on_Z] );
    elseif ismatrix( data )
        inz = double( data(:,1) );
        inaonz = double( data(:,2) );
    end
        
    if nargin == 1
        nsigma = 2;
        zcut = 50.3;
    elseif nargin == 2
        zcut = 50.3; 
    end
    
    
    [hinaonz, binaonz] = hist( inaonz(find( inz < zcut )), ...
                               [min(inaonz):0.001:max(inaonz)] );
    [inaonz_m, inaonz_im] = xb_multigaus_find( hinaonz, 'triglevel', max(hinaonz)/50 );
    pees_init = zeros( 3*numel( inaonz_m ), 1 );
    pees_init(1:3:end) = inaonz_m;
    pees_init(2:3:end) = binaonz( inaonz_im );
    pees_init(3:3:end) = 0.005;
    
    [peesAonZ, ~, errAonZ] = xb_multigaus_fit( [binaonz; hinaonz], pees_init );
    
    isotope = peesAonZ(2:3:end);
    topesigma = peesAonZ(3:3:end);

    peesZ = {};
    errZ = {};
    blob_ellispes = [];
    for ii=1:numel( isotope )
        selection = find( inaonz > isotope(ii) - topesigma(ii) & ...
                          inaonz < isotope(ii) + topesigma(ii) ); 
        if isempty( selection )
            continue;
        end
        [hinz, binz] = hist( inz(selection), [min(inz):0.05:max(inz)] );
        [inz_m, inz_im] = xb_multigaus_find( hinz, 'triglevel', max( hinz )/50 );
        pees_init = zeros( 3*numel( inz_m ), 1 );
        pees_init(1:3:end) = inz_m;
        pees_init(2:3:end) = binz( inz_im );
        pees_init(3:3:end) = 0.2;
        
        [pz, ~, ez] = xb_multigaus_fit( [binz; hinz], pees_init );
        peesZ(end+1) = pz;
        errZ(end+1) = ez;
        
        isotone = pz(2:3:end);
        tonesigma = pz(3:3:end);
        
        for tt=1:numel(isotone)
            A = [ isotone(tt) - nsigma*tonesigma(tt); isotope(ii) + nsigma*topesigma(ii) ];
            B = [ isotone(tt) + nsigma*tonesigma(tt); isotope(ii) + nsigma*topesigma(ii) ];
            C = [ isotone(tt) + nsigma*tonesigma(tt); isotope(ii) - nsigma*topesigma(ii) ];
            D = [ isotone(tt) - nsigma*tonesigma(tt); isotope(ii) - nsigma*topesigma(ii) ];
            blob_ellispes = [blob_ellispes; xb_get_ell_from_rect( [A,B,C,D] )];
        end
    end
end
        
    
    
