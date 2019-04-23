%This function will calculate coulomb cross section, given a beam energy and some cross sections.
%If the GDR/GQR is targeted, the function will know what to do.
%
% [dxs_C/dE, [components]] = xb_xcc( beam_beta, gamma_energy, xs_gam_eX );
%
%parameters:
% -- beam_beta: it's clearly 42.
% -- gamma_energy: some sort of information about from where to where the xcc is calc'd.
% -- xs_gam_eX: an array of functionals representing dxs/dE for different multipolarities.
%returns:
% -- the differential cross section for the given beam energy

function [xcc, components] = xb_xcc( b_beam, e_gam, xs_gam_eX )
	if ischar( xs_gam_eX )
		if strcmp( xs_gam_eX, 'GDR' )
			xs_gam_eX = @__xs_GDR_gam;
			[nfb, ~] = xb_virtualphotons( e_gam*1e3, b_beam, 50, 132, 82, 208 );
		elseif strcmp( xs_gam_eX, 'GQR' )
			xs_gam_eX = @__xs_GQR_gam;
			[~, nfb] = xb_virtualphotons( e_gam*1e3, b_beam, 50, 132, 82, 208 );
		else
			error( 'Unknown string.' );
		end
	end
	
	nb_crux = numel( xs_gam_eX );
	if nb_crux > 2
		error( 'Do not really know how to calculte vPhot number for order > 2...' );
	end
	
	for ii=1:nb_crux
		if ~is_function_handle( xs_gam_eX )
			error( 'Need function handles here' );
		end
	end
	
	e_gam = e_gam(:);
	nfb = nfb(:);
	components = xs_gam_eX( e_gam, b_beam );
	
	xcc = sum( 1./e_gam.*components.*nfb, 2 );
end

function sigma = __xs_GDR_gam( nrg, bt )
	%this notation is not meant for efficiency, but for clarity.
	%I also need to speak MeV, not KeV.
	e_m = @( a ) (31.2*a^(-1/3) + 20.6*a^(-1/6));
	Gamma = @( a ) 0.026*e_m( a )^1.91;
	TRK_e1 = @( z, a ) 60*(a-z)*z/a;
	
	sigma = 2/(pi*Gamma( 132 ))*TRK_e1( 50, 132 )./(1+((nrg.^2-e_m( 132 )^2)./(nrg*Gamma( 132 ))).^2);
	sigma = sigma(:);
end

function [sigma_iv, sigma_is] = __xs_GQR_gam( nrg )
	%here too: efficiency boo, clarity yee
	e_m = @( a ) 64.0*a^-(1/3);
	Gamma = @( a ) 16.0*a^-(1/3);
	TRK_e2is = @( z, a ) 2.2e-4*z^2*a^-(1/3);
	TRK_e2iv = @( z, a ) 2.2e-4*(a-z)*z*a^-(1/3);
	
	sigma_is = 2/(pi*Gamma( 132 ))*TRK_e1is( 50, 132 )*nrg.^2./(1+((nrg.^2-e_m( 132 )^2)./(nrg*Gamma( 132 ))).^2);
	sigma_iv = 2/(pi*Gamma( 132 ))*TRK_e1iv( 50, 132 )*nrg.^2./(1+((nrg.^2-e_m( 132 )^2)./(nrg*Gamma( 132 ))).^2);
	
	sigma_is = sigma_is(:);
	sigma_iv = sigma_iv(:);
end
