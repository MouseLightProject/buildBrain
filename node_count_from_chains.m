function result = node_count_from_chains(chains)
    result = length(unique([chains{:}])) ;
end
