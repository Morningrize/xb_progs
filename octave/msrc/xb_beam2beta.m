%A little utility to convert beta to beam energy (hopefully)
%
% beta = beam2beta( beam_amev, beam_a, beam_z )

function beta = xb_beam2beta( beam_amev, beam_a, beam_z )
	__proton_mass = 938.2720814;
	__neutron_mass = 939.5654135;
	
	rest_m = __proton_mass*beam_z + __neutron_mass*( beam_a - beam_z );
	e = beam_amev*beam_a + rest_m;
	p = sqrt( e.^2 - rest_m^2 );
	beta = p./e;
end
