%this utility converts a reaction probability into a cross section,
%given the target information
%
% xs = xb_rp2xs( rp, target_thickness, target_density_or_material )
%
%arguments:
% --rp : the reaction probabilities to convert. Can be an array (or a spectrum).
% --target_thickness : the target target_thickness, g/cm^2
% --target_density_or_material: either a bivector (density,molar mass
%                               for the density or a string,
%                               hoping the function knows the material.
%NOTE: the density should be given in g/cm^3.

function xs = xb_rp2xs( rp, thk, density_or_material )
    if ischar( density_or_material )
        if strcmp( density_or_material, 'lead' )
            trho = 11.34;
            mmass = 207.2;
        elseif strcmp( density_or_material, 'carbon' )
            trho = 2.267;
            mmass = 12.011;
            warning( 'Using GRAPHITE, assuming nuke physicists too cheap for diamonds.' );
        elseif strcmp( density_or_material, 'graphite' )
            trho = 2.267;
            mmass = 12.011;
        elseif strcmp( density_or_material, 'diamond' )
            trho = 3.515;
            mmass = 12.011;
        elseif strcmp( density_or_material, 'plasitc' )
            trho = 0.926;
            mmass = 2*12.011 + 4*1.008;
            warning( 'Assuminc MDPE (C2H4)_n, density --> 0.926 g/cm^3' );
        else
            error( 'Unkown material requested. Use numbers instead or edit me.' );
        end
    elseif isvector( density_or_material )
        trho = density_or_material(1);
        mmass = density_or_material(2);
    else
        error( 'density_or_material must be either a bivector or a known string.' );
    end
    
    N_a = 6.02214076*1e23; %Brand new Avogadro's constant!
    
    if iscell( rp )
        xs = cell( size( rp ) );
        for ii=1:length( rp )
            xs(ii) = rp{ii}.*(mmass/(thk*N_a)*1e24); %already in barnZ.
        end
    else
        xs = rp.*(mmass/(thk*N_a)*1e24); %already in barnZ.
    end
end
    
    
