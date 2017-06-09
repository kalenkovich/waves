function G = get_gain_matrix(G_bst, channels_idx, ...
                             vertices_idx)
    % From the bst forward model structure extract the forward model
    % matrix, then go from unconstrained to constrained sources, then leave
    % only the required channels and vertices
    % G_bst - bst forward model structure
    G = bst_gain_orient(G_bst.Gain, G_bst.GridOrient);
    G = G(channels_idx, vertices_idx);
    
end