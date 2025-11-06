% Sample data
x = 1:10;
y1 = rand(1,10);
y2 = rand(1,10);
y3 = rand(1,10);

% Plot with different colors (e.g., for variables)
hold on
h1 = plot(x, y1, 'r-');  % red, solid
h2 = plot(x, y2, 'g--'); % green, dashed
h3 = plot(x, y3, 'b:');  % blue, dotted
hd1 = plot(NaN, NaN, 'k-');  % solid
hd2 = plot(NaN, NaN, 'k--'); % dashed
hd3 = plot(NaN, NaN, 'k:');  % dotted
hold off

% First legend: Color legend
legend1 = legend([h1, h2, h3, hd1, hd2,hd3], {'Var A', 'Var B', 'Var C', 'Var D', 'Var F', 'Var G'});
set(legend1, 'Location', 'northeast');

