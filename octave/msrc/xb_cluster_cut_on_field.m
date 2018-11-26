%this function selects cluster based on the content of a field.
%usage: [klz, nb_removed] = xb_cluster_cut_on_nrg( klz, op_handle, field_name )
%       where "op_handle" is a function handle that
%           -takes an aray of energies as argument
%           -returns an array of boolean values (or something that can
%            go into an if statement)
%       NOTE: those for which TRUE is returned are KEPT!

function [klz, nb_removed] = xb_cluster_cut_on_field( klz, op_handle, field_name )
	%klz is a tructure representing the event
	%this will be maintained
	%checkin the input
	if ~is_function_handle( op_handle )
		error( "Second argument **MUST** be a function handle!" );
	end
	if ~ischar( field_name )
		error( 'Third argument **MUST** be a string.' );
	end
	
	if isempty( klz )
		nb_removed = 0;
		return;
	end
	
	%do the parallel execution
    proc_args = cell(2);
    proc_args(1) = op_handle;
    proc_args(2) = field_name;
	[klz, nb_removed] = xb_parproc( klz, @_processor, proc_args );
end

%the usual processor function
function [klz, nb_removed] = _processor( klz, proc_args )
    op_handle = proc_args{1};
    field_name = proc_args{2};

	%loop clear the uninsteresting stuff.
	nb_removed = sum( [klz.n] );
	for ii=1:length( klz )
		if isfield( klz(ii).clusters, field_name )
			keep_idx = find( op_handle( [klz(ii).clusters.( field_name )](:) ) );
			klz(ii).clusters = klz(ii).clusters( keep_idx );
		
			%at the end of things, multiplicity update
			klz(ii).n = length( klz(ii).clusters );
		else
			warning( ['At index', num2str(ii), ...
			          ' the requested field "', ... 
			          field_name, ' wasn`t found.'] );
		end
	end
	nb_removed = nb_removed - sum( [klz.n] );
end
