%% load independence polynomial for fully-filled lattice

data_SG = load_full_lattice('../data/independence_polynomial_square_grid/independence_polynomial_square_lattice', (10:30)');
data_K = load_full_lattice('../data/independence_polynomial_kings_graphs/independence_polynomial_diagonal_coupled', (10:28)');



%% functions
function data = load_full_lattice(file_string, L_ind)
data = struct;
for ind = 1:numel(L_ind)
    L = L_ind(ind);
    fid = fopen(sprintf([file_string, '_n%d.dat'],L));
    ind_poly = sym(textscan(fid, '%s'));
    fclose(fid);
    
    no_ind_sets = sum(ind_poly);
    entropy_constant = double(no_ind_sets^(1/L^2));
    
    data(ind).ind_poly = ind_poly;
    data(ind).L = L_ind(ind);
    data(ind).no_ind_sets = no_ind_sets;
    data(ind).entropy_constant = entropy_constant;
end
end

























