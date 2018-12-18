%this function should be able to do all the parallel setup and lounching of
%cutomized processors.
%
% [dataset_proc, nb_removed] = xb_parproc( dataset, processor, proc_args )
%
%Parameters:
% -- dataset: the dataset (an XB kind) on which to operate.
% -- processor: the function that will actually do the work on dataset
% -- proc_args: a cell array with all the proper arguments to the processor funciton
%               which clearly should accept a cell array as argument and be able
%               to parse it.
%Returns:
% -- dataset_proc: the processed dataset
% -- nb_removed: the number of removed elements

function [dataset, nb_removed] = xb_parproc( dataset, processor, proc_args )
    %prepare to split the events in the number of processes
    dataset = dataset(:)';
	nb_events = numel( dataset );
	adj = mod( nb_events, nproc );
	if adj adj = nproc - adj; end
	idx_part = 1:(nb_events+adj)/nproc:(nb_events+adj);
	idx_part = linspace( idx_part, ...
	                     idx_part+(nb_events+adj)/nproc-1, ...
	                     (nb_events+adj)/nproc );
	
	%split the events into a cell
	dataset_part = {};
	dataset_rest = [];
	nbr_rest = 0;
	nb_proc = 1;
	for ii=1:nproc
		try
			dataset_part(ii) = dataset(idx_part(ii,:));
		catch
            %do padding
            dataset_part(ii) = [dataset(idx_part(ii):end),dataset(1:adj)];
		end
		nb_proc += 1;
	end

	%do the parallel execution
	proc_handle = @( p ) processor( p, proc_args );
	try
		pkg load parallel;
		[dataset_part, nb_removed_part] = parcellfun( nb_proc, proc_handle, ...
		                                              dataset_part, 'VerboseLevel', 0 );
	catch
		warning( 'Parallel package not available. This will take a while.' );
		[dataset_part, nb_removed_part] = cellfun( proc_handle, dataset_part, ...
                                                   'UniformOutput', false );
        nb_removed_part = cell2mat( nb_removed_part );
	end
	
	%stitch together the stuff
	dataset = reshape( dataset_part, [], 1 );
    if iscell( dataset ) dataset = cell2mat( dataset ); end
	nb_removed = sum( nb_removed_part ); 

	%remove the padding and prune the empty ones
	dataset = dataset(1:nb_events);
	dataset = dataset( find( [dataset.n] ) );
end
