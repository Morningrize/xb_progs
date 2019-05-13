%this function calculates a photoabsorption cross section from an oscillator strength
%"strength" as defined B( El, J_i->J_f ) = 1/(2Ji-1)*|<J_f||M(El)||J_i>|^2
%
% phxs = xb_strength2phxb( energy, fsd_nrg, strength, l );
%
%arguments:
% energy : either an array or a value, the excitation energy to target. KeV.
% fsd_nrg : the final state density as a function of the excitation energy.
%           Either an array with the same dimensions of energy or a function handle.
% strength : the strength in question. The unit of measure expected is e^2*b^l
% l : the multipolarity order.
%NOTE: the kind (magnetic, electric) of the strength --> the returned cross section
%      is not an explicit parameter: it comes with the number of the strength.
%returns:
% phxs: the photoabsorption cross section corresponding to the given energy(es)
%       and strength
%NOTE: this function does not sum over final states, but considers one single B(pl, i->f )
%      at one time. If you need to sum, eval more time and then sum for phxs^pl.

function phxs = xb_strength2phxc( nrg, fsd, strength, l )
    if is_function_handle( fsd )
        fsd = fsd( nrg );
    elseif size( fsd ) ~= size( nrg )
        error( 'fsd_nrg and energy MUST have the same dimension.' );
    end
    
    _hbar = 1.05457180013e-34; %reduced Planck constant, J*s
    _c = 2.99792458e23; %speed of light, fm/s
    
    nrg = nrg*1.6021766208e-16; %KeV to J
    
    strength = strength*100^l; %from b^l to fm^(2*l) 
    
    ki = nrg/(_hbar*_c); %this thing, as it is, has units of 1/fm.
    phxs = (2*pi)^3*(l+1)/(l*xb_semifact( 2*l+1 )^2) * ...
           fsd.*ki.^(2*(l-1))*strength;
    %And this has the units of fm^2. Which is correct.
    %NOTE: it stinks SO MUCH like there have been a e^-2 left off in Ber87 in this formula
    %      (Dominic had a parenthesis off, ki^(2*l-1) <-- ki^(2*(l-1)) ).
    
    phxs /= 100; %fm^2 to barns.
    
end
    
    
