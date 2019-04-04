%This function returns a distribution df/dE for the virtual photon number exchanged
%By two nuclei flying past each other.
%
% [dN/dE_E1, dN/dE_E2] = xb_virtualphotons( energy, proj_velocity, proj_mass, target_charge, target_mass )
%
%parameters:
% -- energy: it's the energy OF THE PHOTONS EMITTED, KeV.
% -- proj_velocity: it's the beta of the incoming ion.
% -- target_charge: it's the Z number of the charge.
% -- proj_mass: projectile's mass number
% -- target_charge: target's Z number
% -- target_mass: target's mass number
%returns:
% -- a matrix, representing the number of virtual photons exchanged
%    at the given conditions for E1 and E2

function [nbf_e1, nbf_e2] = xb_virtualphotons( nrg, bt, Ap, Zt, At )
    if ~isscalar( bt ) && ~isscalar( nrg )
        error( 'Cannot handle array of velocities and gamma energies, yet.' );
    end

    _alpha = 0.0072973525664; %the fine structure constant, CODATA14
    _hbar = 6.62607015e-34/(2*pi); %reduced Plank constant.
    _c = 299792458; %speed of light
    _gamma = 1./sqrt( 1 - bt^2 ); %lorentz factor
    _b_min = 1.34*( Ap^(1/3) + At^(1/3) -0.75*( Ap^(-1/3) + At^(-1/3) ) )*1e-15; %minimum impact parameter, SI units.
    _bees = 1.75e5*rand( size( nrg ) )*1e-15; %randomly generated impact parameters (m) up to atomic radius (more than that should mean closer to another nucleus, in a very simplistic fashion).
    _bees( find( _bees < _b_min ) ) = 0; %this will set _xi = 0 --> nbf = 0 for that entry.
    _E_max = _hbar*_c*bt*_gamma/(_b_min*1.6021766208e-16) %maximum excitation energy (closer is a crash). This comes out in KeV;
    
    if ~isempty( find( nrg > _E_max ) )
		warning( 'Some energies are not accessible. Setting them to 0' );
		nrg( find( nrg > _E_max ) ) = 0; %This will also null the adiabaticity parameter
	end
    
    _xi = (bt.^-1)*nrg*_b_min/(_hbar*_c*_gamma ); %the adiabaticity parameter, assuming all collisions at _b_min ATM.
    _K0 = besselk( 0, _xi );
    _K1 = besselk( 1, _xi );
    _norm = 2*Zt^2*_alpha./(pi*bt^2);

    
    nbf_e1 = _norm*( _xi.*_K0.*_K1 - _xi.^2*bt^2/2.*( _K1.^2 - _K0.^2 ) );
    
    nbf_e2 = _norm*( 2*(1 - bt.^2).*_K1.^2 + _xi.*(2 - bt.^2).*_K0.*_K1 ...
            - _xi.^2*bt.^4/2.*(_K1.^2 - _K0.^2));
end
