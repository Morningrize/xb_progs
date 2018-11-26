%this function selects cluster based on their total energy.
%usage: [klz, nb_removed] = xb_cluster_cut_on_nrg( klz, op_handle )
%       where "op_handle" is a function handle that
%           -takes an aray of energies as argument
%           -returns an array of boolean values (or something that can
%            go into an if statement)
%       NOTE: those for which TRUE is returned are KEPT!

function [klz, nb_removed] = xb_cluster_cut_on_nrg( klz, op_handle )
	%klz is a tructure representing the event
	%this will be maintained
	if ~is_function_handle( op_handle )
		error( "Second argument **MUST** be a function handle!" );
	end
	if isempty( klz )
		nb_removed = 0;
		return;
	end
	
	proc_args = cell(1); proc_args(1) = op_handle;
	[klz, nb_removed] = xb_parproc( klz, @_processor, proc_args );
end

%the processor function
function [klz, nb_removed] = _processor( klz, proc_args )
    op_handle = proc_args{1};
    
	%loop clear the uninsteresting stuff.
	nb_removed = sum( [klz.n] );
	for ii=1:length( klz )
		keep_idx = find( op_handle( [klz(ii).clusters.sum_e](:) ) );
		klz(ii).clusters = klz(ii).clusters( keep_idx );
		
		%at the end of things, multiplicity update
		klz(ii).n = length( klz(ii).clusters );
	end
	nb_removed = nb_removed - sum( [klz.n] );
end
