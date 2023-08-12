function diff_mat = coord2diff(node_vals)
   
for node_i = 1:size(node_vals,1)
    for node_j = 1:size(node_vals,1)
    diff_mat(node_i,node_j) = 0.5 * (node_vals(node_i,1) - node_vals(node_j,1)).^2;
    end
end
end