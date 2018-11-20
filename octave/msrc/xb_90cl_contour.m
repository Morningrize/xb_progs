%this tiny function returns the centroid and the 90% cl contour of the data
%assuming they distribute according to a 2D gaussian.
%
% ell = xb_3sigma_contour( data_x, data_y )
%
%NOTE: the multiplier for one sigma has been found on
%      https://people.richland.edu/james/lecture/m170/tbl-chi.html
%      Better ones are welcome.
%parameters:
% -- datax, data_y: x and y of your data
%returns:
% -- ell: the three sigma ellipse contour. It's a structure containing:
%         ell.ctr -- the centroid
%         ell.a   -- x// semiaxis
%         ell.b   -- y// semiaxis
%         ell.rot -- the rotation of the ellipse
%

function ell = xb_90cl_contour( data_x, data_y )
	%assuming a normal distribution, mean and covrinace should be as good as a fit
	%if not better.
	ell.ctr = mean( [data_x(:), data_y(:)] );
	c_m = cov( [data_x(:), data_y(:)] );
	
	%now, the contour
	sigma_mult = sqrt( 4.605 );
	[e_vec, e_val] = eig( c_m ); %eigenvalues and vectors
	ell.a = sigma_mult*sqrt( e_val(1,1) ); %ellipse's A (for 3sigma)
	ell.b = sigma_mult*sqrt( e_val(2,2) ); %ellipse's B (for 3sigma)
	ell.rot = atan2( e_vec(2,1), e_vec(1,1) ); %ellipse's rotation
end
