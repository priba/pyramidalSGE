function h = draw_graphs_clustering_from_adjmat(A, cluster)

N = size(A,1);

vertices = [cos(2*pi*(1:N)/N); sin(2*pi*(1:N)/N)]';

classes = unique(cluster);
k = length(classes);

colors=hsv(k);               % plots points with appropriate colors
colormap(colors)


[x,y] = gplot(A, vertices);
h = figure; hold on;
plot(x, y, '-', 'LineWidth', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'MarkerSize', 15 );
for i = 1:length(classes)
    idx = classes(i)==cluster;
    plot(vertices(idx,1), vertices(idx,2), 'ok', 'LineWidth', 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', colors(i,:), 'MarkerSize', 15 );
end;
hold off;
axis off;

end