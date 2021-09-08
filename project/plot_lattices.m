% Plot some lattices 

rng('default');
L = 10; % lattice dimension
p = 0.8; % filling factor
[x,y] = ind2sub([L,L],1:L^2);
xy_lattice_full = [x',y'];

index_site_remove = randperm(L^2);
index_site_remove = index_site_remove(1:round((1-p)*L^2));
xy_lattice_partial = xy_lattice_full;
xy_lattice_partial(index_site_remove,:) = [];

%% plots 
markersize = 40;
linewidth = 2;

%% Square lattice, full
[edges, ~, ~] = UnitDiskGraph(xy_lattice_full); % radius = 1

figure; hold on; axis equal; axis tight; axis off;
G = graph(edges(:,1), edges(:,2));
h = plot(G, 'XData', xy_lattice_full(:,1), 'YData', xy_lattice_full(:,2), 'LineWidth',linewidth, 'MarkerSize',markersize, 'EdgeColor','k', 'NodeColor', 'k');
h.NodeLabel = {};
% saveas(h,'../paper/figures/square_lattice.eps');

%% Square lattice, partial
[edges, ~, ~] = UnitDiskGraph(xy_lattice_partial); % radius = 1

figure; hold on; axis equal; axis tight; axis off;
G = graph(edges(:,1), edges(:,2));
h = plot(G, 'XData', xy_lattice_partial(:,1), 'YData', xy_lattice_partial(:,2), 'LineWidth',linewidth, 'MarkerSize',markersize, 'EdgeColor','k', 'NodeColor', 'k');
h.NodeLabel = {};
% saveas(h,'../paper/figures/square_lattice_partial.eps');


%% Square diagonal lattice, full
[edges, ~, ~] = UnitDiskGraph(xy_lattice_full/(sqrt(2)+0.01)); % radius = sqrt(2)

figure; hold on; axis equal; axis tight; axis off;
G = graph(edges(:,1), edges(:,2));
h = plot(G, 'XData', xy_lattice_full(:,1), 'YData', xy_lattice_full(:,2), 'LineWidth',linewidth, 'MarkerSize',markersize, 'EdgeColor','k', 'NodeColor', 'k');
h.NodeLabel = {};
% saveas(h,'../paper/figures/square_diagonal.eps');


%% Square diagonal lattice, partial
[edges, ~, ~] = UnitDiskGraph(xy_lattice_partial/(sqrt(2)+0.01)); % radius = sqrt(2)

figure; hold on; axis equal; axis tight; axis off;
G = graph(edges(:,1), edges(:,2));
h = plot(G, 'XData', xy_lattice_partial(:,1), 'YData', xy_lattice_partial(:,2), 'LineWidth',linewidth, 'MarkerSize',markersize, 'EdgeColor','k', 'NodeColor', 'k');
h.NodeLabel = {};
% saveas(h,'../paper/figures/square_diagonal_partial.eps');


