function ForFigure = getFigData(numA, ForFigure)
if ~exist('ForFigure', 'var')
    tmp = gcf;
    ForFigure(numA).dat = tmp.UserData.currTsData;
    ForFigure(numA).dat.veText = tmp.UserData.ui.r2Text.String;  
%     ForFigure(numA).dat.veText = tmp.UserData.ui.r2Text.String;  
end
    tmp = gcf;
    ForFigure(numA).dat = tmp.UserData.currTsData;
    ForFigure(numA).dat.veText = tmp.UserData.ui.r2Text.String;     
end