%this function returns (in principle) the virtual photon number for an arbitrary multipolarity
%and magnetic/electric mode.
%
% dN/dE = xb_VPN( kind, order, speces )
%
%parameters:
% -- kind: magnetic 'M' or electric 'E' transition.
% -- order: order of the transition above
% -- specs: a structure containint  (at lease) these fields:
%    -- gam_e: the energy (sampled range, in an array) of the emitted gammas
%    -- bt: the beta
%    -- Zp : the charge of the projectile
%    -- Ap : the mass number of the projectile
%    -- Zt : the charge of the target
%    -- At : the mass number of the target
%    NOTE: any additional field will be ignored.
%returns:
% -- the distribution of these photons (not normalized!) in terms of the energy.

function vpn = xb_VPN( kind, order, specs )
	if ~isscalar( bt ) && ~isscalar( nrg )
        error( 'Cannot handle array of velocities and gamma energies, yet.' );
    end


    _alpha = 0.0072973525664; %the fine structure constant, CODATA14
    _hbar = 6.62607015e-34/(2*pi); %reduced Plank constant.
    _c = 299792458; %speed of light, m/s
    _gamma = 1./sqrt( 1 - specs.bt^2 ); %lorentz factor
    _b_min = 1.34*( specs.Ap^(1/3) + specs.At^(1/3) ...
    -0.75*( specs.Ap^(-1/3) + specs.At^(-1/3) ) )*1e-15; %minimum impact parameter, SI units.
    _bees = 1.75e5*rand( size( specs.nrg ) )*1e-15; %randomly generated impact parameters (m) up to atomic radius (more than that should mean closer to another nucleus, in a very simplistic fashion).
    _bees( find( _bees < _b_min ) ) = 0; %this will set _xi = 0 --> nbf = 0 for that entry.
    _E_max = _hbar*_c*specs.bt*_gamma/(_b_min); %maximum excitation energy (closer is a crash). This is in Joules now;
    
    nrg = nrg*1.6021766208e-16; %KeV to J
    if ~isempty( find( nrg > _E_max ) )
		warning( 'Some energies are not accessible. Setting them to 0' );
		nrg( find( nrg > _E_max ) ) = 0; %This will also null the adiabaticity parameter
    end
    
    _xi = (specs. bt.^-1)*nrg*_b_min/(_hbar*_c*_gamma ); %the adiabaticity parameter.
    
    l = order;
    m = [-l:l];
    vpn = specs.Zp^2*_alpha*
end

%support function, contains bessel functions
function val = __g( m, xi )
	K_m = besselk( m, xi );
	K_mp1 = besselk( m+1, xi );
	
	
	val = pi*xi.^2*( K_mp1.^2 - K_m.^2 - 2*m./xi.*K_m.*K_mp1 );
end

%
function val = __G_mag( m, l, x )
	if x >= 1
		warning( '__G_el defined only for x < 1.' );
		val = NaN;
		return;
	end
	
	imp = find( m >= 0 );
	imn = find( m < 0 );
	if ~isempty( imn )
		val = (-1)^abs( m(imn)+1 )*__G_mag( abs( m(imn) ), l, x );
	end
	
	m = m(imp);
	val = [val, i^(l+m+1)*sqrt( 16*pi )/factorial( factorial( l*(2*l + 1 ) ) ) ...
	*sqrt( factorial( l-m )/factorial( l+m ) )*sqrt( x.^2 -1 ) ...
	*m*legendre( m, l, x )];
end

function val = __G_el( m, l, x )
	if x >= 1
		warning( '__G_el defined only for x < 1.' );
		val = NaN;
		return;
	end
	
    imp = find( m >= 0 );
	imn = find( m < 0 );
	if ~isempty( imn )
		val = (-1)^abs( m(imn)+1 )*__G_el( abs( m(imn) ), l, x );
	end	
	%note that these are vectors (if x is a scalar, as it should).
	m = m(imp);
	lef_lm1 = legendre( m, l-1, x );
	lef_lp1 = legendre( m, l+1, x );
	
	val = i^(l+m)*sqrt( 16*pi )/factorial( factorial( l*(2*l +1) ) ) ...
	*sqrt( factorial( l-m )./factorial( l+m ) )*sqrt( x.^2 -1 ) ...
	.*( (l+1)*(l+m)/(2*l +1).*lef_lm1 - l*(l-m+1)/(2*l+1).*lef_lp1 );
end
