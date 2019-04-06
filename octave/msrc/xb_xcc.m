%This function will calculate coulomb cross section, given a beam energy and some cross sections.
%If the GDR/GQR is targeted, the function will know what to do.
%
% dxs_C/dE = xb_xcc( beam_energy, gamma_energy, xs_gam_eX );
%
%parameters:
% -- beam_energy: it's clearly 42.
% -- gamma_energy: some sort of information about from where to where the xcc is calc'd.
% -- xs_gam_eX: an array of functionals representing dxs/dE for different multipolarities.
%returns:
% -- the differential cross section for the given beam energy

function xcc = xb_xcc( e_beam, e_gam, xs_gam_eX )
	if ischar( xs_gam_eX )
		if strcmp( xs_gam_eX, 'GDR' )
			xs_gam_ex = @__xs_GDR_gam;
		elseif strcmp( xs_gam_eX, 'GQR' )
			xs_gam_eX = @__xs_GQR_gam;
		else
			error( 'Unknown string.' );
		end
	end
	
	nb_crux = numel( xs_gam_eX );
	if nb_crux > 2
		error( 'Do not really know how to calculte vPhot number for order > 2...' );
	end
	
	for ii=1:nb_crux
		if ~is_function_handle( xs_gam_eX(ii) )
			error( 'Need function handles here' );
		end
	end
	
	xcc = %...
end

function __xs_GDR_gam( nrg )
	%...
end

function __xs_GQR_gam( nrg )
	%...
end
