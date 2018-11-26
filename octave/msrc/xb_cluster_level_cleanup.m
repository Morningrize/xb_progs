%this function implement Andrea Horvat's minimum reasonable level cut
%It operates on clusters
%
% [klz_clean, nb_removed] = xb_cluster_level_cleanup( klz, level_energy )
%
% parameters:
% --klz: a cluster struct array
% --level_energy: the target level energy.

function [klz, nb_removed] = xb_cluster_level_cleanup( klz, level_energy )
	if isempty( klz )
		nb_removed = 0;
		return;
	end
	
	proc_args = cell(1);
    proc_args(1) = level_energy;
    [klz, nb_removed] = xb_parproc( klz, @_processor, proc_args );
    
end

function [klz, nb_removed] = _processor( klz, level_energy )
    nb_removed = sum( [klz.n] );
    for ii=1:numel( klz )
        if ~isempty( find( [klz(ii).clusters.sum_e] > level_energy ) )
            continue;
        else
            klz(ii).n = 0;
            klz(ii).clusters = [];
        end
    end
    nb_removed = nb_removed - sum( [klz.n] );
end
