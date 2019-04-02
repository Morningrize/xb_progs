%This function returns a distribution df/dE for the virtual photon number exchanged
%By two nuclei flying past each other.
%
% dN/dE = xb_virtualphotons_E1( energy, proj_velocity, proj_mass, target_charge, target_mass )
%
%parameters:
% -- energy: it's the energy
% -- proj_velocity: it's the beta of the incoming ion.
% -- target_charge: it's the Z number of the charge.
%returns:
% -- a number, representing the number of virtual photons exchanged
%    at the given conditions.

function nbf = xb_virtualphotons_E1( nrg, b, Ap, Zt, At )
    if isvector( nrg ) && isvector( b ) && ( numel( nrg ) ~= numel( b ) )
        error( 'If energy and beta are vectors, they must have the same size' );
    end

    _alpha = 0.0072973525664; %the fine structure constant, CODATA14
    _hbar = 6.62607015e-34/(2*pi); %reduced Plank constant.
    _c = 299792458; %speed of light
    _gamma = 1/.sqrt( 1 - b^2 ); %lorents factor
    _b_min = 1.34*( Ap^(1/3) + At^(1/3) -0.75*( Ap^(-1/3) + At^(-1/3) ) ); %minimum impact parameter.
    _E_max = _hbar*_c*b*_gamma/_b_min; %maximum excitation energy (closer is a crash)
    _xi = nrg*_b_min./(_h_bar*_c*b*_gamma ); %the adiabaticity parameter, assuming all collisions at _b_min ATM.
    _K0 = besselk( 0, _xi );
    _K1 = besselk( 1, _xi );
    
    if( E > _E_max )
        nbf = 0;
        return;
    end
    
    nbf = 2*Zt^2*_alpha/*(pi*b^2)*( _xi.*_K0.*_K1 - _xi.^2*b^2/2.*( _K1.^2 - _K0.^2 ) );
end
