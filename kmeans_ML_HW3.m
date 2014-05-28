clc;
clear all;
close all;
startTime=tic;

data1 = csvread('data1.csv');
data2 = csvread('data2.csv');
data3 = csvread('data3.csv');

figure(1);
plot(data1(:,1),data1(:,2),'ob');
% figure(2);
% plot(data2(:,1),data2(:,2),'b');
% figure(3);
% plot(data3(:,1),data3(:,2),'b');

% number of clusters
K = 5;


numberOfTrainingData = size(data1,1);
numberOfFeatureDimensions = size(data1,2);

for index=1:200,


% initialize cluster centers
randCenterIndices = randi([1 numberOfTrainingData], 1, K);
clusterCenters = data1(randCenterIndices,:);

data1Label = zeros(numberOfTrainingData, 1);
data1DistanceFromClusterCenter = Inf(numberOfTrainingData, K);

% Assign each data point to a cluster: label each data point with a cluster
% number

changeInLabelsFlag = 1;

while changeInLabelsFlag == 1,
    
    changeInLabelsFlag = 0;
    
    for i = 1:numberOfTrainingData,
        for j = 1:K,
           data1DistanceFromClusterCenter(i,j) = EuclideanDistance(data1(i,:),clusterCenters(j,:));
        end
        % check to see if any data point's label is updated
        [~, newLabel] = min(data1DistanceFromClusterCenter(i,:));
        if  data1Label(i) ~= newLabel,
            changeInLabelsFlag = 1;
        end
        % update label maybe?
        [~, data1Label(i)] = min(data1DistanceFromClusterCenter(i,:));
        
    end
    
    
    % Re-estimate cluster centers
    
    Mu = Inf(K, numberOfFeatureDimensions);
    
    for i = 1:K,
        
        Mu(i,:) = mean(data1(find(data1Label == i),:));
        
        minDiff = EuclideanDistance(data1(1,:), Mu(i,:));
        
        for j = 2:numberOfTrainingData
            if minDiff > EuclideanDistance(data1(j,:), Mu(i,:)),
                minDiff = EuclideanDistance(data1(j,:), Mu(i,:));
                closestDataPointToMuIndex = j;
            end
        end
        
        clusterCenters(i,:) = data1(closestDataPointToMuIndex,:);
    end
    
end


hold on;
plot(clusterCenters(:,1),clusterCenters(:,2),'oc');
hold on;
    
    
end

title({'K-Means clustering with K = 5, for 200 times'});
legend('Original Data Points', 'Cluster Centers');
xlabel('X');
ylabel('Y');

% save as PNG and EPS
saveas(1, 'fig_Kmeans_with_random_initialization', 'epsc2');
saveas(1, 'fig_Kmeans_with_random_initialization', 'png');
saveas(1, 'fig_Kmeans_with_random_initialization', 'fig');
saveas(1, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Kmeans_with_random_initialization'], 'epsc2');

% save the environment variables
save ([datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS')]);

runningTime = toc(startTime);
runningTime


