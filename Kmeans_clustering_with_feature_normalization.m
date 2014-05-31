% K-means clustering algorithm
% Authors: Amir Azarbakht and Mandana Hamidi
% azarbaam@eecs.oregonstate.edu
% 2014-05-30

clc;
clear all;
close all;
startTime=tic;

data3 = csvread('data3.csv');

% number of clusters
K = 2;

numberOfTrainingData = size(data3,1);
numberOfFeatureDimensions = size(data3,2);

% Normalize the data set
for counterFeatures=1: numberOfFeatureDimensions  
  data3Normalized(:,counterFeatures)= (data3(:,counterFeatures)-mean(data3(:,counterFeatures)))/std(data3(:,counterFeatures)) ;
end

figure(1);
plot(data3Normalized(:,1),data3Normalized(:,2),'ob');

% for the numberOfRuns, here = 200
for index=1:200,
    
    
    % initialize cluster centers
    
%     randCenterIndices = randi([1 numberOfTrainingData], 1, K);
    % use randperm to get *unique* random numbers
    randCenterIndices = randperm(numberOfTrainingData, K);
    
    
    clusterCenters = data3Normalized(randCenterIndices,:);
    
    data3NormalizedLabel = zeros(numberOfTrainingData, 1);
    data3NormalizedDistanceFromClusterCenter = Inf(numberOfTrainingData, K);
    
    % Assign each data point to a cluster: label each data point with a cluster
    % number
    
    changeInLabelsFlag = 1;
    
    while changeInLabelsFlag == 1,
        
        changeInLabelsFlag = 0;
        
        for i = 1:numberOfTrainingData,
            for j = 1:K,
                data3NormalizedDistanceFromClusterCenter(i,j) = EuclideanDistance(data3Normalized(i,:),clusterCenters(j,:));
            end
            % check to see if any data point's label is updated
            [~, newLabel] = min(data3NormalizedDistanceFromClusterCenter(i,:));
            if  data3NormalizedLabel(i) ~= newLabel,
                changeInLabelsFlag = 1;
            end
            % update label maybe?
            [~, data3NormalizedLabel(i)] = min(data3NormalizedDistanceFromClusterCenter(i,:));
            
        end
        
        
        % Re-estimate cluster centers
        
        Mu = Inf(K, numberOfFeatureDimensions);
        
        for i = 1:K,
            
            Mu(i,:) = mean(data3Normalized(find(data3NormalizedLabel == i),:));
            
            
            minDiff = EuclideanDistance(data3Normalized(1,:), Mu(i,:));
            
            for j = 2:numberOfTrainingData
                if minDiff > EuclideanDistance(data3Normalized(j,:), Mu(i,:)),
                    minDiff = EuclideanDistance(data3Normalized(j,:), Mu(i,:));
                    closestDataPointToMuIndex = j;
                end
            end
            
%             clusterCenters(i,:) = data3Normalized(closestDataPointToMuIndex,:);
            clusterCenters(i,:) = Mu(i,:);
            % Part 1. b.
            
            
            
        end
        
    end
    
       
    
    hold on;
    plot(clusterCenters(:,1),clusterCenters(:,2),'*c');
    hold on;
    
    
    withinClusterSum = zeros(K,1);
    
    % clustering done
    
    % calculate withing cluster sum of distances
    for i = 1:K,
        
        %      clusterI = zeros(size((data3Normalized(find(data3NormalizedLabel == i),:))));
        clusterI = data3Normalized(find(data3NormalizedLabel == i),:);
        for j = 1:size(clusterI,1),
            withinClusterSum(i) = withinClusterSum(i) + (EuclideanDistance(Mu(i,:), clusterI(j,:)))^2;
        end
        
    end
    
    if min(withinClusterSum) == 0
       test = 0; 
    end
    
    MinMaxMeanSD(index, 1) = min(withinClusterSum); 
    MinMaxMeanSD(index, 2) = max(withinClusterSum);
    MinMaxMeanSD(index, 3) = mean(withinClusterSum);
    MinMaxMeanSD(index, 4) = std(withinClusterSum);
    
    
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
saveas(1, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Kmeans_with_random_initialization'], 'fig');
saveas(1, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Kmeans_with_random_initialization'], 'png');

% plot the metrics
figure(2);
hold on;
plot((1:1:200), MinMaxMeanSD(:,1), '+g');
hold on;
plot((1:1:200), MinMaxMeanSD(:,2), '*b');
hold on;
% to plot the standard deviation(s), we use errorbar on the mean(s)
% errorbar((1:5:200), MinMaxMeanSD(1:5:end,3), MinMaxMeanSD(1:5:end,4), '.r');
errorbar((1:1:200), MinMaxMeanSD(1:1:end,3), MinMaxMeanSD(1:1:end,4), '.r');

hold on;
plot((1:1:200), MinMaxMeanSD(:,3), 'ok');
hold on;
title({'Min, Max, Mean and STD of Within-Cluster Sum of Squared Distances for K-Means clustering with K = 5, for 200 times'});
legend('Min', 'Max', 'Standard deviation', 'Mean');
hold on;
xlabel('200 runs');
ylabel('Within-Cluster Sum of Squared Distances');

saveas(2, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Within_cluster_sum_of_squared_distances_200runs'], 'epsc2');
saveas(2, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Within_cluster_sum_of_squared_distances_200runs'], 'fig');
saveas(2, ['fig_' datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS') '_' 'Within_cluster_sum_of_squared_distances_200runs'], 'png');

% save the environment variables
save ([datestr(date, 'YYYY-mm-dd') '_' datestr(now, 'HH-MM-SS')]);

runningTime = toc(startTime);
% print run time
runningTime


