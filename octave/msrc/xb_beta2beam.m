%A little utility to convert beta to beam energy (hopefully)
%
% beam_amev = beta2beam( beta, beam_a, beam_z )

function beam_amev = xb_beta2beam( beta, beam_a, beam_z )
	__proton_mass = 938.2720814;
	__neutron_mass = 939.5654135;
	
	rest_m = __proton_mass*beam_z + __neutron_mass*( beam_a - beam_z );
	gamma = sqrt( 1 - beta.^2 ).^-1;
	beam_amev = ( gamma - 1 ).*rest_m/beam_a;
end
