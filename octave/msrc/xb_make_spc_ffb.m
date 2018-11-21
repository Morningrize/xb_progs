%this routine produces a comparison between full, front and back spectrum
%of some data set ({a,}data or cluster) for the various interesting section of
%the CB
%
% [hst_ffb, binz_ffb, herr_ffb, nrg_ffb] = xb_make_spc_ffb( dataset, bin, [figure, options] )
%
%parameters:
% --dataset: the dataset (data or clusters) to act upon.
% --bin: either the bin width or an array representing the binning.
% --options:
% ----a figure: a figure object to draw upon
% ----'singles': it's the default and does a singles spectrum
% ----'sum': produces a sum spectrum instead.
% ----'field': followed by another string, do a spectrum of that field.
%
%returns:
% --hst_ffb: a cell array with the full, front and back spectrum.
% --binz_ffb: the bin arrays for those specrta
% --herr_ffb: the errors on the counts per bin
% --nrg_ffb: the array content of the field processed.

function [hst_ffb, binz_ffb, herr_ffb, nrg_ffb] = xb_make_spc_ffb( dataset, bin, varargin )
	if isfield( dataset, 'clusters' )
		iswhat = 'cluster';
	elseif isfield( dataset, 'e' ) && isfield( dataset, 'he' )
		iswhat = 'data';
	else
		error( 'dataset is neither cluster nor data' );
	end
	
	spiccer_data = @xb_data_nrg;
    spiccer_klz = @xb_cluster_nrg;
	
	if nargin == 2 && nargout == 0
		fig = figure;
	else
        for aa=1:numel( varargin )
            if isfigure( varargin{aa} )
                fig = varargin{aa};
            elseif ischar( varargin{aa} ) && strcmp( varargin{aa}, 'singles' )
                spiccer_data = @xb_data_nrg;
                spiccer_klz = @xb_cluster_nrg;
            elseif ischar( varargin{aa} ) && strcmp( varargin{aa}, 'sums' )
                spiccer_data = @xb_data_sumnrg;
                spiccer_klz = @xb_cluster_sumnrg;
            elseif ischar( varargin{aa} ) && strcmp( varargin{aa}, 'field' )
                spiccer_data = @( p ) xb_data_field( p, varargin{aa+1} );
                spiccer_klz = @( p ) xb_cluster_field( p, varargin{aa+1} );
            end
	end
	
	hst_ffb = cell( 3, 1 );
	binz_ffb = cell( 3, 1 );
	herr_ffb = cell( 3, 1 );
	nrg_ffb = cell( 3, 1 );
	icbf = [xb_ball_neigh( 81, 5 ).i];
	[~, icbb] = setdiff( [xb_ball_at( [1:162] ).i], icbf );
	ohf = @(p) xb_op_cbi( p, icbf );
	ohb = @(p) xb_op_cbi( p, icbb );
	
	if strcmp( iswhat, 'data' )
		front = xb_data_cut_on_field( dataset, ohf, 'i' );
		back = xb_data_cut_on_field( dataset, ohb, 'i' );
		nrg_ffb(1) = spiccer_data( dataset );
		nrg_ffb(2) = spiccer_data( front );
		nrg_ffb(3) = spiccer_data( back );
	elseif strcmp( iswhat, 'cluster' )
		front = xb_cluster_cut_on_field( dataset, ohf, 'centroid_id' );
		back = xb_cluster_cut_on_field( dataset, ohb, 'centroid_id' );
		nrg_ffb(1) = spiccer_klz( dataset );
		nrg_ffb(2) = spiccer_klz( front );
		nrg_ffb(3) = spiccer_klz( back );
	end
	clear dataset front back;
	
	[hst_ffb(1), binz_ffb(1), herr_ffb(1)] = xb_make_spc( nrg_ffb{1}, bin );
	[hst_ffb(2), binz_ffb(2), herr_ffb(2)] = xb_make_spc( nrg_ffb{2}, bin );
	[hst_ffb(3), binz_ffb(3), herr_ffb(3)] = xb_make_spc( nrg_ffb{3}, bin );
	
	if exist( 'fig' ) && isfigure( fig );
		figure( fig );
		stairs( binz_ffb{1}, hst_ffb{1}, 'linewidth', 2 ); hold on;
		stairs( binz_ffb{2}, hst_ffb{2}, 'linewidth', 2 );
		stairs( binz_ffb{3}, hst_ffb{3}, 'linewidth', 2 );
		set( gca, 'fontsize', 24, 'linewidth', 2, 'yscale', 'log' );
		grid on;
		if isscalar( bin ); bstep = bin;
		else bstep = bin(2)-bin(1); end
		ylabel( ['#/',num2str(bstep),'KeV'] );
		xlabel( 'KeV' );
		legend( { 'full CB', 'front CB', 'back CB' } );
	end
end
