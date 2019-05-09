%this function calculates a photoabsorption cross section from an oscillator strength
%"strength" as defined B( El, J_i->J_f ) = 1/(2Ji-1)*|<J_f||M(El)||J_i>|^2
%
% phxs = xb_strength2phxb( energy, fsd_nrg, strength, l );
%
%arguments:
% energy : either an array or a value, the excitation energy to target. KeV.
% fsd_nrg : the final state density as a function of the excitation energy.
%           Either an array with the same dimensions of energy or a function handle.
% strength : the strength in question.
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
    
    _ech = -1.602176620899*1e-19; %charge of the electron, C
    _hbar = 6.62607015e-34/(2*pi); %reduced Planck constant, J*s
    _c = 299792458; %speed of light, m/s
    
    nrg = nrg*1.6021766208e-16; %KeV to J
    
    ki = nrg/(_hbar*_c);
    phxs = (2*pi)^3*(l+1)/(l*xb_semifact( 2*l+1 )^2) * ...
           fsd.*ki.^(2*l-1)*strength*_ech^2;
    %and it plops out in barns. Hopefully.
    
end
    
    
