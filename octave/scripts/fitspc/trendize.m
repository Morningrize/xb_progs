%A tiny function to plot the trend of a model (without fmincon)
%
% [minrg, minval] = trendize( span, model, fig )

function [minrg, minval] = trendize( span, model, fig )
	if nargin == 2
		fig = figure;
	else
		figure( fig );
		hold on;
	end

	tr = [];
	for ee=span(1):span(2)
		tr = [tr, model( ee )];
	end
	
	plot( [span(1):span(2)], tr, 'x', 'linewidth', 2 );
	set( gca, 'fontsize', 24, 'linewidth', 2 );

	[minval, ii] = min( tr );
	minrg = [span(1):span(2)](ii); grid on;
end
