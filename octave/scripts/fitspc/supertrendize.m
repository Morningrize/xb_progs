%A tiny function to plot the trend of a bivariate model
%
% tr = supertrendize( span, model, fig )

function tr = supertrendize( span, model, nsteps, fig )
	if nargin == 2
        fig = figure;
        nsteps = 100;
	elseif nargin == 3
		fig = figure;
	else
		figure( fig );
		hold on;
	end

	stepping = [ (span(2)-span(1))/nsteps, (span(4)-span(3))/nsteps ];
	x = span(1):stepping(1):span(2);
	y = span(3):stepping(2):span(4);
	tr = zeros( length( x ), length( y ) );
	for ii=1:length( x )
        for jj=1:length( y )
            tr(ii,jj) = model( [x(ii), y(jj)] );
        end
	end
	
	plot3( x, y, tr, 'rx', 'linewidth', 2 );
	set( gca, 'fontsize', 24, 'linewidth', 2 );
    grid on;
end
