%this function selects hits based on their energy.
%usage: [evt, nb_removed] = xb_data_cut_on_nrg( evt, op_handle )
%       where "op_handle" is a function handle that
%           -takes an aray of energies as argument
%           -returns an array of boolean values (or something that can
%            go into an if statement)
%       NOTE: those for which TRUE is returned are KEPT!

function [evt, nb_removed] = xb_data_cut_on_nrg( evt, op_handle )
	if ~is_function_handle( op_handle )
		error( "Second argument **MUST** be a function handle!" );
	end
	if isempty( evt )
		nb_removed = 0;
		return;
	end
	
	proc_args = cell(1); proc_args(1) = op_handle;
	[evt, nb_removed] = xb_parproc( evt, @_processor, proc_args );
end

%processor function
function [evt, nb_removed] = _processor( evt, proc_args )
    op_handle = proc_args{1};
    
	%loop clear the uninsteresting stuff.
	nb_removed = sum( [evt.n] );
	for ii=1:length( evt )
		keep_idx = find( op_handle( [evt(ii).e](:) ) );
		if size( evt(ii).i ) evt(ii).i = evt(ii).i( keep_idx ); end
		if size( evt(ii).e ) evt(ii).e = evt(ii).e( keep_idx ); end
		if size( evt(ii).he ) evt(ii).he = evt(ii).he( keep_idx ); end
		if size( evt(ii).t ) evt(ii).t = evt(ii).t( keep_idx ); end
		if size( evt(ii).pt ) evt(ii).pt = evt(ii).pt( keep_idx ); end
		
		%at the end of things, multiplicity update
		evt(ii).n = length( keep_idx );
		evt(ii).sum_e = sum( [evt(ii).e] );
	end

	nb_removed = nb_removed - sum( [evt.n] );
end
