clear;
clc;

filename =('Lumbar_FlexExt_Kinematics.xlsx');%load kinematic file
Matrix = xlsread(filename);%create time-series data as a matrix.
%Data is nXm matrix (n = time series; m = participant)

%k = the number of clusters
%change k variable to change the number of clusters calculated
k = 10;

%Generates a random number within the size of the matrix for each cluster
%These random numbers will be used as the initial centroids for each
%cluster. It is possible that the random number generator will select the
%same number twice. When this happens, just run it again.
Cluster_Rand = zeros(1,k);
for ii = 1:k
    Cluster_Rand(1,ii) = randi(size(Matrix,2));
end

%Here I find the time-series that matches the random numbers
%These time-series will be used as the initial centroids to start off
%k-means analysis
Cluster_Centroid = zeros(length(Matrix),k);
for ii = 1:k
    Cluster_Centroid(:,ii) = Matrix(:,Cluster_Rand(ii));
end

%This while loop will run k-means clustering until no time-series is reassigned to a different cluster.
%when the "Cluster_Num_Difference" is 0 (no difference in clusters between previous and current interations of k-means)
%k-means analysis will stop. 
Cluster_Num_Difference_Ave = 1;%I chose one for this variable as an arbitrary number just to create the variable before the while loop.
n = 0; %This variable will tell me how many time k-means runs before it stops.
while Cluster_Num_Difference_Ave ~= 0
for ii = 1:k
    Cluster_DTW(ii,:) = zeros(1,size(Matrix,2));
end

%determines the distance (using dtw) between the initial centroids and all the vectors in the matrix 
for ii = 1:k
    for a=1:size(Matrix,2)
    Cluster_DTW(ii,a) = dtw(Cluster_Centroid(:,ii),Matrix(:,a));
    end
end

%in order for me to group all the Cluster_DTW variables together, I need to round the data.
for ii = 1:k
   Cluster_DTW(ii,:) = round(Cluster_DTW(ii,:));
end

%Here I find which centroid each vector is closest to by finding the
%minimum DTW distance from all three centroids
Cluster_min = double(1:size(Matrix,2));
for ii = 1:size(Matrix,2)
    Cluster_min(ii) = min(Cluster_DTW(:,ii));
end

%Now I use the minimum distance found above and will match that distance
%measure to the appropriate vector
for ii = 1:k
    Cluster_Num(ii,:) = zeros(1,size(Matrix,2));
    Cluster_length(ii) = length(find(Cluster_DTW(ii,:)==Cluster_min));
    Cluster_Num(ii,1:Cluster_length(1,ii)) = find(Cluster_DTW(ii,:)==Cluster_min);
    Cluster_Num_1{ii} = {Cluster_Num(ii,:)};
end

%removes zeros from cell arrays
for ii = 1:k
    Cluster_Num_1{1,ii}{1,1}(Cluster_Num_1{1,ii}{1,1} == 0) = [];
end

%Here I will take those find the vectors for each cluster within the Matrix
for ii = 1:k
    Cluster_Vectors_Doub = Matrix(:,Cluster_Num_1{1,ii}{1,:});
    Cluster_Vectors_Cell{ii} = {Cluster_Vectors_Doub};
    Cluster_Vectors_Doub = [];
end

%Sum all the vectors within each Cluster 
for ii = 1:k
    Cluster_Sum{ii} = 0;
    for a = 1:length(Cluster_Num_1{1,ii}{1,1})
        Cluster_Sum{ii} = Cluster_Sum{ii} + Cluster_Vectors_Cell{1,ii}{1,1}(:,a);
    end
end

%Average each cluster
for ii = 1:k
    Cluster_Ave{ii} = Cluster_Sum{ii}/length(Cluster_Num_1{1,ii}{1,1});
end

%Find how many participants are in each cluster
for ii = 1:k
    Cluster_length(ii) = length(Cluster_Num_1{1,ii}{1,1}(1,:));
end

%This will compare the subject numbers from the previous run of k-means to
%the current run of k-means to determine if the clusters changed. If the
%subjects within each cluster do not change, k-means will stop.
%This for loop will only start after the first complete cycle of the while
%loop.
if n > 0
    for ii = 1:k
    Cluster_Num_Ave(ii) = mean(Cluster_Num_1{1,ii}{1,1});
    Cluster_Num_Previous_Ave(ii) = mean(Cluster_Num_Previous{1,ii}{1,1});
    Cluster_Num_Difference(ii) = Cluster_Num_Ave(ii) - Cluster_Num_Previous_Ave(ii);
    end
    Cluster_Num_Difference_Ave = mean(Cluster_Num_Difference);
end

%change variable with subject numbers to compare to the next run of k-means
Cluster_Num_Previous = Cluster_Num_1; 

%I now change Cluster average variables to cluster centroid variables in
%order to replace the intial centroids at the beginning of the while loop
for ii = 1:k
    Cluster_Centroid(:,ii) = Cluster_Ave{1,ii}(:,1);
end
n = n+1;
end

%This will be used to display the number of participants in each cluster on
%the graph 
for ii = 1:k
    Cluster_Total(1:2,ii) = [ii Cluster_length(1,ii)];
end

%Graph the final results
figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
hold on
for ii = 1:k
    plot(Cluster_Centroid(:,ii),'LineWidth',3.0)
end
hold on
xline(42,'--');
xline(64,'--');
xline(80,'--');
title('K-means++ Clustering Trunk Flexion and Extension', 'FontSize', 16)
set(get(gca,'title'),'Position',[50 370 1.00011])
xlabel('Trunk Flexion and Extension Cycle (%)','FontSize',13,'FontWeight','bold')
set(get(gca,'xlabel'),'Position',[50 -110 1.00011])
ylabel('Trunk Flexion/Extension Angle (°)','FontSize',13,'FontWeight','bold')
hold on
xlim([0 100])
annotation('textbox',[.11 .075 .2 .01],'String','Neutral','EdgeColor','none', 'FontSize', 12)
annotation('textbox',[.42 .075 .2 .01],'String','Max Flexion','EdgeColor','none', 'FontSize', 12)
annotation('textbox',[.60 .075 .2 .01],'String','Neutral','EdgeColor','none', 'FontSize', 12)
annotation('textbox',[.71 .075 .2 .01],'String','Hyperextension','EdgeColor','none', 'FontSize', 12)
annotation('textbox',[.885 .075 .2 .01],'String','Neutral','EdgeColor','none', 'FontSize', 12)
%legend(sprintfc('Cluster %d (n = %d)',Cluster_Total),'Location','east', 'FontSize', 12)