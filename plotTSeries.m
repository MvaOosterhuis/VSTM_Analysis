
% load('ForFigure.mat')
colvec = {'r', 'b', 'g'};
pvec = {'.r', '.b', '.g'};
figure('Position', [50, 50, 1500, 600]); 
plot(ForFigure(1).dat.x, zeros(size(ForFigure(1).dat.x)), 'k', 'LineWidth', 2);
hold on;
axis([0 106.5 -30 30])
colrbar = zeros(size(ForFigure(1).dat.x, 3));
for n = 1:6
    colrbar(1+6*(n-1):6+6*(n-1),:) = 10+10*(n-1);
    colrbar(37+6*(n-1):42+6*(n-1),:) = 60-10*(n-1);
end
image([0 106.5],[30 -30],colrbar')
% ylim = [-30 25];
hold on; 

for x = 1:3
    plot(ForFigure(x).dat.x, ForFigure(x).dat.pred, colvec{x}, 'LineWidth', 3)
    plot(ForFigure(x).dat.x, ForFigure(x).dat.raw, 'o', 'MarkerFaceColor', colvec{x}, 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'LineWidth', 2)
%     plot(xvec, newdat, 'o', 'MarkerFaceColor', colvec{x}, 'MarkerSize', 8, 'MarkerEdgeColor', 'k')
end
% 
% plot(ForFigure(1).dat.x,zeros(size(ForFigure(1).dat.x)), 'k')
% 
% plot(ForFigure(1).dat.x,ForFigure(1).dat.pred, 'LineWidth', 3); 
% 
% plot(ForFigure(1).dat.x,ForFigure(1).dat.raw, '.', 'MarkerSize', 30);
% axis([0 106.5 -30 25])
% 
% % figure('Position', [50, 50, 1500, 600]); 
% plot(ForFigure(2).dat.x,ForFigure(2).dat.pred, 'LineWidth', 3); 
% hold on; 
% plot(ForFigure(2).dat.x,ForFigure(2).dat.raw, '.', 'MarkerSize', 30);
% axis([0 106.5 -30 25])
% 
% plot(ForFigure(3).dat.x,ForFigure(3).dat.pred, 'LineWidth', 3); 
% hold on; 
% plot(ForFigure(3).dat.x,ForFigure(3).dat.raw, '.', 'MarkerSize', 30);
% axis([0 106.5 -30 25])

%%

newsize = length(ForFigure(1).dat.x)/2;
xvec = ForFigure(1).dat.x(1:newsize);
figure('Position', [50, 50, 800, 600]);
hold on;
plot(xvec,zeros(size(xvec)), 'k', 'LineWidth', 2.5)
axis([0 53.5 -30 25])
colvec = {'r', 'b', 'g'};
pvec = {'.r', '.b', '.g'};
colrbar = zeros(size(xvec, 3));
for n = 1:6
    colrbar(1+6*(n-1):6+6*(n-1),:) = 10+10*(n-1);
%     colrbar(6:12,:) = 2;
%     colrbar(31:36,:) = 3;
end
image([0 53.5],[-30 30],colrbar')

for x = 1:3
    tmp = ForFigure(x).dat.pred(1:newsize);
    tmp2 = flip(ForFigure(x).dat.pred(newsize+1:end));
    newdat = (tmp+tmp2)/2;
    plot(xvec, newdat, colvec{x}, 'LineWidth', 3)
    
    tmp = ForFigure(x).dat.raw(1:newsize);
    tmp2 = flip(ForFigure(x).dat.raw(newsize+1:end));
    newdat = (tmp+tmp2)/2;
    plot(xvec, newdat, 'o', 'MarkerFaceColor', colvec{x}, 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'LineWidth', 2)
%     plot(xvec, newdat, 'o', 'MarkerFaceColor', colvec{x}, 'MarkerSize', 8, 'MarkerEdgeColor', 'k')
end
