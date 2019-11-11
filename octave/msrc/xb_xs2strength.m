%This function is the reverse (hopefully) of xb_strength2phxs.
%
% strength = xb_xs2strength( nrg, fsd, xs, l )
%
%arguments:
% nrg : the excitation energy(ies) at which this is happening
% fsd : the final state density as a function of excitation energy.
%       Just put one in it...
% xs : it's the cross section. It's the photoabsorbtion one, not the coulex one.
% l : the multipolarity order you're interested in.
%returns:
% strength : the strength corresponding to the cross section passed. 
%            It can be a vector, but only one multiploarity (shouldn't be a problem)

function strength = xb_xs2strength( nrg, fsd, xs, l )
    if is_function_handle( fsd )
        fsd = fsd( nrg );
    elseif size( fsd ) ~= size( nrg )
        error( 'fsd and erg must have the same dimension.' );
    end

    _hbar = 1.05457180013e-34; %reduced Planck constant, J*s
    _c = 2.99792458e23; %speed of light, fm/s

    nrg = nrg*1.6021766208e-16; %KeV to J
    xs *= 100; %barns to fm^2
    ki = nrg/(_hbar*_c); %this thing, as it is, has units of 1/fm.

    strength = ( (2*pi)^3*(l+1)/(l*xb_semifact( 2*l+1 )^2) * ...
               fsd.*ki.^(2*(l-1)) ).^-1.*xs;

    strength /= 100^l; %fm^(2*l) to b^l
end
