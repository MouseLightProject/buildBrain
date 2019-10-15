function assignment_matrix = neuron_branch_assignment(neurons,branches)
    % graph: the complete graph structure
    % assignment_matrix - a sparse matrix 
    n_neurons = length(neurons);
    n_branches = length(branches);
    proximity_threshold = 3;

    branch_2_neurons = zeros(n_branches, n_neurons);

    for neuron_index = 1:n_neurons
        disp(neuron_index)
        neuron = neurons(neuron_index);
        if isfield(neuron, 'GT') && ~isempty(neuron.GT)
            neuron_pix_locs = neuron.GT.swcpixlocs;
            parfor branch_index = 1:n_branches
                branch_pix_locs = branches(branch_index).subs;
                [branch_length, ~] = size(branch_pix_locs);
                branch_2_neuron_dist = pdist2(neuron_pix_locs, branch_pix_locs, 'euclidean');
                [~, bpi] = find(branch_2_neuron_dist < proximity_threshold);
                ubpi = unique(bpi);
                matched_length = length(ubpi);
                if matched_length > 0.7 * branch_length
                    % only keep it if the branch matches the GT more than
                    % 70%
                    branch_2_neurons(branch_index, neuron_index) = length(ubpi);
                end
            end
        end
    end
    
    assignment_matrix = sparse(branch_2_neurons);
end