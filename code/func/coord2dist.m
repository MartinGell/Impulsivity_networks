function dist_mat = coord2dist(node_coords)
   
for node_i = 1:size(node_coords,1)
    dist_mat(node_i,:) = sqrt(sum(bsxfun(@minus, node_coords(node_i,:), node_coords) .^ 2, 2));
end
end