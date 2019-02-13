numOfSubs = 19;
front = {};
parietal = {};
temporal = {};
data = {};
timepoints = 116;
T = {};
%T = timepoints * ones(numOfSubs,1);
cd '/Users/shivamdwivedi/Desktop//control/control_20'
for i = 1:numOfSubs
    temp = importdata(i + ".txt");
    normTemp = transpose(temp);
    data{i,1} = temp(1:timepoints,:);
    T{i,1} = timepoints;
end
% for i = 1:numOfSubs
%     temp = importdata(i + ".txt");
%     normTemp = transpose(temp);
%     data{i,1} = normTemp.data(1:timepoints,:);
%     T{i,1} = timepoints;
% end

% for i = 1:numOfSubs
%     temp = importdata(i + ".txt");
%     normTemp = transpose(temp);
%     front{i,1} = normTemp.data(1:timepoints,3:16);
%     parietal{i,1} = normTemp.data(1:timepoints,59:62);
%     temporal{i,1} = normTemp.data(1:timepoints,81:90);
%     T{i,1} = timepoints;
% end