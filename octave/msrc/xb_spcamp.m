%this utility takes in a spectrum and "amplifies it" --i.E. it will treat it as a probability distribution
%and produce a random number generator with it, and then produce a spectrum with more (or less) counts in
%the integral, as requested. Smoothing of the original spectrum may happen.
%
% [aspc, aspcerr] = xb_spcamp( spc, bin, amp_factor, [smooth] )
%
%parameters:
% -- spc: the original spectrum
% -- bin: the binning of the spectrum
% -- amp_factor: the factor with which to multiply the original integral ( < 1 is allowed )
% -- smooth: optionally, an integer, odd factor for a Savitzky-Golay filetrin of the original spectrum
%returns:
% -- aspc: the amplified spectrum
% -- spcerr: the error on the amplified spectrum (relative error is conserved)

function [aspc, aspcerr] = xb_spcamp( spc, bin, amp_factor, smooth )
    intg = sum( spc );
    aintg = round( amp_factor*intg ); %has to be an integer, of course
    if nargin == 4
        try
            pkg load signal
        catch
            error( 'I need the package "signal" to do smoothing' );
        end
        smooth = round( smooth );
        if mod( smooth, 2 ) == 0; smooth += 1; end
        spc = sgolayfilt( spc, 2, smooth );
        intg = sum( spc );
    end
    erel = spc.^(-1/2);

    spc /= intg;
    gen = xb_arbitrary_distro( spc, bin );

    deposits = gen( 1, aintg );
    aspc = hist( deposits, bin );
    aspcerr = aspc.*erel;
end
