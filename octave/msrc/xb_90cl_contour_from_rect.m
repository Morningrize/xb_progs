%this tiny function returns the centroid and the 90% cl contour of the data
%inside a rectangle, assuming they distribute according to a 2D gaussian.
%
% [ell, idx_in] = xb_90cl_contour_from_rect( data_x, data_y, rect )
%
%NOTE: the multiplier for 3sigma has been found on
%      https://people.richland.edu/james/lecture/m170/tbl-chi.html
%      Better ones are welcome.
%parameters:
% -- datax, data_y: x and y of your data
% -- rect:  an array describing a rectangle (4 scalars for an unrotated one, 2x4 else)
%           NOTE: if it's not a rectangle, all of YOUR life has no meaning
%
%The rectancle can be rotated, but the ordering of the vertices must be as follows:
%   A B C D
%  [x x x x;
%   y y y y]
%
%  C-------B
%  |       |
%  |       |
%  D-------A
%returns:
% -- ell: the three sigma ellipse contour. It's a structure containing:
%         ell.ctr -- the centroid
%         ell.a   -- x// semiaxis
%         ell.b   -- y// semiaxis
%         ell.rot -- the rotation of the ellipse
%

function [ell, idx_in] = xb_90cl_contour_from_rect( datax, datay, rect )
    big_ell = xb_get_ell_from_rect( rect );
    idx_big = xb_is_in_ellipse( datax, datay, big_ell );
    ell = xb_90cl_contour( datax(idx_big), datay(idx_big) );
    if nargout == 2
        idx_in = xb_is_in_ellipse( datax, datay, ell );
    end
end
    
