%this function implement Andrea Horvat's minimum reasonable level cut
%It operates on clusters
%
% [klz_clean, nb_removed] = xb_cluster_level_cleanup( klz, level_energy )
% [klz_clean, nb_removed] = xb_cluster_level_cleanup( klz, level_energy, cmp_sign );
%
% parameters:
% --klz: a cluster struct array
% --level_energy: the target level energy.
% --cmp_sign: the sign of comparison: '>' or '<'

function [klz, nb_removed] = xb_cluster_level_cleanup( klz, level_energy, cmp_sign )
	if isempty( klz )
		nb_removed = 0;
		return;
	end
	
	proc_args = cell(1);
    proc_args(1) = level_energy;
    if nargin == 2 || nargin == 3 && cmp_sign == '>'
        [klz, nb_removed] = xb_parproc( klz, @_processor_greater, proc_args );
    elseif nargin == 3 && cmp_sign == '<'
        [klz, nb_removed] = xb_parproc( klz, @_processor_lesser, proc_args );
    else
        error( 'Invalid comparison sign' );
    end
    
end

%that's not elegant, but I can't be bothered at the moment.
function [klz, nb_removed] = _processor_greater( klz, proc_args )
    level_energy = proc_args{1};
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

function [klz, nb_removed] = _processor_lesser( klz, proc_args )
    level_energy = proc_args{1};
    nb_removed = sum( [klz.n] );
    for ii=1:numel( klz )
        if isempty( find( [klz(ii).clusters.sum_e] > level_energy ) )
            continue;
        else
            klz(ii).n = 0;
            klz(ii).clusters = [];
        end
    end
    nb_removed = nb_removed - sum( [klz.n] );
end
