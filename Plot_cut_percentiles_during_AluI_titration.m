load('Cut_ratios_all_sites_111P_4.mat', 'fcut_all_sites', 'c_vector')

Quantiles = quantile(fcut_all_sites', [0.05:0.05:0.95]);
noQuantiles = size(Quantiles, 1);
for q = 1:noQuantiles
    Quantiles(q,:) = transpose(smooth(Quantiles(q,:), 21));
end

x = c_vector';
C = [1, .6, .6];
dC = ([1, .95, .95] - C) / 8;

figure
clear p
hold all
for q = 1:2:5
    p(round((q+1)/2)) = patch([x, fliplr(x)],[Quantiles(20-q,:), fliplr(Quantiles(q,:))], C+(9-q)*dC, 'EdgeColor', 'none');
end
l = plot(x, Quantiles(10,:), 'r', 'linewidth', 2);

ylim([0 1])
xlim([0 max(c_vector)])

xlabel('AluI concentration (nM)')
ylabel('Cut fraction')

h = legend(p, {'5%-95%', '15%-85%', '25%-75%'}, 'location', 'EO');
v = get(h,'title');
set(v,'string','Percentiles');

set(gca, 'YGrid', 'on', 'GridLineStyle', '--', 'layer', 'top')

set(gcf, 'PaperPositionMode', 'auto')
print(gcf, '-depsc', '-painters', 'Cut_ratio_percentiles.eps');
