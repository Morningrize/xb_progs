%this function selects hits based on the content of a field.
%usage: [evt, nb_removed] = xb_trk_cut_on_field( evt, op_handle )
%       where "op_handle" is a function handle that
%           -takes an aray of energies as argument
%           -returns an array of boolean values (or something that can
%            go into an if statement)
%       NOTE: those for which TRUE is returned are KEPT!

function [evt, nb_removed] = xb_trk_cut_on_field( evt, op_handle, field_name )
	if ~is_function_handle( op_handle )
		error( "Second argument **MUST** be a function handle!" );
	end
	if ~ischar( field_name )
		error( 'Third argument **MUST** be a string.' );
	end
	if ~isfield( evt, field_name )
		warning( ['No field named "',field_name,'" in the given structure.'] );
		nb_removed = 0;
		return;
	end
	if isempty( evt )
		nb_removed = 0;
		return;
	end
	%do t
	he parallel execution
    proc_args = cell(2);
    proc_args(1) = op_handle;
    proc_args(2) = field_name;
	[evt, nb_removed] = xb_parproc( evt, @_processor, proc_args );
end

%processor function: does the actaul cutting.
%it's called in parallel... for a bit of speed.
function [evt, nb_removed] = _processor( evt, proc_args )
    op_handle = proc_args{1};
    field_name = proc_args{2};
    
	%loop clear the uninsteresting stuff.
	nb_removed = sum( [evt.n] );
	for ii=1:length( evt )
		keep_idx = find( op_handle( [evt(ii).( field_name )](:) ) );
		if size( evt(ii).fragment_A ) evt(ii).fragment_A = evt(ii).fragment_A( keep_idx ); end
		if size( evt(ii).fragment_Z ) evt(ii).fragment_Z = evt(ii).fragment_Z( keep_idx ); end
		if size( evt(ii).incoming ) evt(ii).incoming = evt(ii).incoming( keep_idx ); end
		if size( evt(ii).outgoing ) evt(ii).outgoing = evt(ii).outgoing( keep_idx ); end
		if size( evt(ii).fragment_beta )
			evt(ii).fragment_beta = evt(ii).fragment_beta( keep_idx );
		end

		%at the end of things, multiplicity update
		evt(ii).n = length( keep_idx );
	end

	nb_removed = nb_removed - sum( [evt.n] );
end
