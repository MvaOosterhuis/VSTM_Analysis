
colvec = {'r', 'b', 'g'};
pvec = {'.r', '.b', '.g'};
figure('Position', [50, 50, 1500, 600]); 
plot(ForFigure(1).dat.x, zeros(size(ForFigure(1).dat.x)), 'k', 'LineWidth', 2);
hold on;
axis([0 106.5 -1 1])
colrbar = zeros(size(ForFigure(1).dat.x, 3));
for n = 1:6
    colrbar(1+6*(n-1):6+6*(n-1),:) = 10+10*(n-1);
    colrbar(37+6*(n-1):42+6*(n-1),:) = 60-10*(n-1);
end
image([0 106.5],[30 -30],colrbar')
% ylim = [-30 25];
% hold on; 
% % figure; hold on;
% for x = 1:3
%     plot(Ms(x).M.currTsData.x, Ms(x).M.currTsData.raw./max(Ms(x).M.currTsData.raw), 'o', 'MarkerFaceColor', colvec{x}, 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'LineWidth', 2)
%     plot(Ms(x).M.currTsData.x, Ms(x).M.currTsData.pred./max(Ms(x).M.currTsData.raw), colvec{x}, 'LineWidth', 3)
%     
% end

% figure;
hold on; 
% figure; hold on;
for x = 1:3
    div = max(max(Ms(x).M.currTsData.raw)-min(Ms(x).M.currTsData.raw),max(Ms(x).M.currTsData.pred)-min(Ms(x).M.currTsData.pred));
    tmp = (Ms(x).M.currTsData.raw)./div;
    tmp2 = (Ms(x).M.currTsData.pred)./div;
    plot(Ms(x).M.currTsData.x, tmp2, colvec{x}, 'LineWidth', 3)
end

for x = 1:3
    div = max(max(Ms(x).M.currTsData.raw)-min(Ms(x).M.currTsData.raw),max(Ms(x).M.currTsData.pred)-min(Ms(x).M.currTsData.pred));
    tmp = (Ms(x).M.currTsData.raw)./div;
    tmp2 = (Ms(x).M.currTsData.pred)./div;
    plot(Ms(x).M.currTsData.x, tmp, 'o', 'MarkerFaceColor', colvec{x}, 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'LineWidth', 2)
end

figure; hold on;
for x = 1:3
[tmp, tmp2] = calculateGaussian(Ms(x).M.rfParams(1), Ms(x).M.rfParams(3), log([0.01:0.01:6]));
tmp=tmp./max(tmp);
plot(exp(tmp2), tmp); axis([0 6 0 1]); axis square
end


%%
% Ms(3).M = M;

%%
% figure; hold on;
% for x = 1:3
% [tmp, tmp2] = calculateGaussian(log(4.95), 1, log([0.01:0.01:6]));
% tmp=tmp./max(tmp);
% plot(exp(tmp2), tmp); axis([0 6 0 1]); axis square
% end