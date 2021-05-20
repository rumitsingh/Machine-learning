%% Silhouette Coefficient
% calculate a = average distance of i to the points in the same cluster
% calculate b = min(average distance of i to points in another cluster)
% Silhouette coefficient = 1 - a/b 
% if a < b, (or s = b/a - 1 if a > or equal to b, not the usual case)
% typically between -1 and 1
% the closer to 1 the better

k=z;%In the Hierarchical Clustering code, the variable for number of clusters is "z"
%Switch z variable to k to work with the code.

%This calculates (a) variable in the silhouette equation.
%This calculates the distance of each vector to other vectors within the
%same cluster.
for ii = 1:k
    n = 1;
    while n<size(Hierarchical_Vectors{1,ii},2)+1
        for a = 1:size(Hierarchical_Vectors{1,ii}(1,:),2)
            Cluster_Vector(ii).Silhouette(a,n) = dtw(Hierarchical_Vectors{1,ii}(:,n),Hierarchical_Vectors{1,ii}(:,a));
        end
        n = n+1;
    end
end

%This sorts the values in each row from smallest to greatest number.
%This is so I can delete the column with all the zeros. 
%Zeros represent the distance between the same vector, which we do not
%need.
for ii = 1:k
    for a = 1:size(Cluster_Vector(ii).Silhouette,2)
        Cluster_Vector(ii).Silhouette(a,:) = sort(Cluster_Vector(ii).Silhouette(a,:));
    end
end

%This deletes the column with the zeros
for ii = 1:k
    Cluster_Vector(ii).Silhouette(:,1) = [];
end

%This calculates the average of for each row vector
for ii = 1:k
    for a = 1:length(Cluster_Vector(ii).Silhouette)
        Cluster_Vector(ii).Silhouette_Ave(a,:) = mean(Cluster_Vector(ii).Silhouette(a,:));
    end
end

%This calculates the distance between each vector within each cluster to
%the centroid of the other clusters. This will help determine the nearest
%centroid each vector is closest two (other than the centroid they are already in)
for ii = 1:k
    for b = 1:k
        for a = 1:size(Hierarchical_Vectors{1,ii}(1,:),2)
            Cluster_Dist(ii).Between(b,a) = dtw(Hierarchical_Vectors{1,ii}(:,a),Matrix(:,b));
        end
    end
end

%The above section also calculated the distance between each vector and
%their own centroid which we do not need. This section will delete those
%rows because we do not need them. 
for ii = 1:k
    Cluster_Dist(ii).Between(ii,:) = [];
end

%This calculates the minimum distance values within each vector distance
%measures.
if k > 2 %this section does not work when the k = 2
    for ii = 1:k
        Cluster_Dist(ii).Min = min(Cluster_Dist(ii).Between);
    end
end

%This section subtracts the minimum values just calculated with the
%Cluster_Dist.Between variable. A zero represents a vector that is closest
%to that centroid.
if k > 2 %this section does not work when the k = 2
    for ii = 1:k
        for a = 1:size(Cluster_Dist(ii).Between,1)
            Cluster_Dist(ii).Sim(a,:) = Cluster_Dist(ii).Between(a,:) - Cluster_Dist(ii).Min;
        end
    end
end

%In order to have this code calculate silhouette coefficients for any
%number of centroids (k), I created the Clus variable, which creates a series
%of numbers based on how many centroid there are.
for ii = 1:k
    Clus(ii,1:k) = 1:k;
end

%This will eliminate the centroid number of that centroid. 
%So within the first column which is centroid 1, the 1 will be changed to
%a zero and centroid two the two will be changed to a zero and so on.
for ii = 1:k
    for a = 1:k
        if Clus(ii,a) == ii
            Clus(ii,a) = 0;
        end
    end
end

%I want to delete the zeros within this variable so I start by sorting the
%rows
for ii = 1:k
    Clus(ii,:) = sort(Clus(ii,:));
end

%This deletes the column with the zeros
Clus(:,1) = [];

%This section will change the Cluster_Dist.Sim variable to show which cluster
%number that vector is closest to. 
if k > 2
    for ii = 1:k
        for a = 1:size(Cluster_Dist(ii).Sim,1) 
            for b = 1:size(Cluster_Dist(ii).Sim(a,:),2)
                if Cluster_Dist(ii).Sim(a,b) == 0
                    Cluster_Dist(ii).Sim(a,b) = Clus(ii,1);
                elseif Cluster_Dist(ii).Sim(a,b) > 0
                    Cluster_Dist(ii).Sim(a,b) = Clus(ii,2);
                end
            end
        end
    end
end

%This section calculates the distance between each individual vector in one
%cluster and to each vector in the other clusters. 
for ii = 1:k
    n = 1;
    while n<size(Hierarchical_Vectors{1,ii},2)+1
        for b = 1:k
            for a = 1:size(Hierarchical_Vectors{1,b},2)
                Silhouette(b,ii).Cluster(n,a) = dtw(Hierarchical_Vectors{1,ii}(:,n),Hierarchical_Vectors{1,b}(:,a));
            end
        end
        n = n+1;
    end
end

%Since I only care about the distances in vectors to the opposite cluster,
%I do not need the calculates of the distances for the same vector. This
%section deletes those sections not needed.
for ii = 1:k
    Silhouette(ii,ii).Cluster = [];
end

%This averages the distance measure calculated for Silhouette.Cluster
%variable when k = 2
if k == 2
    for ii = 1:k
        for a = 1:size(Hierarchical_Vectors{1,ii},2)
            Silhouette_Between(ii).Ave(a,1) = mean(Silhouette(Clus(ii,1),ii).Cluster(a,:));
        end
    end
end

%This averages the distance measure calculated for Silhouette.Cluster
%variable when k is greater than 2
if k > 2
    for ii = 1:k
        for a = 1:size(Hierarchical_Vectors{1,ii},2)
            if Cluster_Dist(ii).Sim(1,a) == Clus(ii,1)
                Silhouette_Between(ii).Ave(a,1) = mean(Silhouette(Clus(ii,1),ii).Cluster(a,:));
            elseif Cluster_Dist(ii).Sim(1,a) == Clus(ii,2)
                Silhouette_Between(ii).Ave(a,1) = mean(Silhouette(Clus(ii,2),ii).Cluster(a,:));
            end
        end
    end
end

%This calculates the silhouette coefficient. This code will calculate the
%equation depending on if a>b or if b>a.
for ii = 1:k
    for a = 1:size(Hierarchical_Vectors{1,ii},2)
        if Cluster_Vector(ii).Silhouette_Ave(a,1)<Silhouette_Between(ii).Ave(a,1)
            Cluster_Vector(ii).SC(a,1) = 1 - (Cluster_Vector(ii).Silhouette_Ave(a,1)/Silhouette_Between(ii).Ave(a,1));
        elseif Cluster_Vector(ii).Silhouette_Ave(a,1)>Silhouette_Between(ii).Ave(a,1)
            Cluster_Vector(ii).SC(a,1) = (Silhouette_Between(ii).Ave(a,1)/Cluster_Vector(ii).Silhouette_Ave(a,1))-1;
        end
    end
end

%This will organize all the silhouette coefficients for each time-series in
%one column when the number of clusters is greater than 2.
if k >2
    Silhouette_Total = [];
for ii = 1:k
    for a = 1:length(Cluster_Vector(ii).SC)
    Silhouette_Total(end+1) = Cluster_Vector(1,ii).SC(a,1);
    end
end
end

%This will organize all the silhouette coefficients for each time-series in
%one column when the number of clusters equal 2.
if k == 2
    Silhouette_Total = [];
for ii = 1:k
    for a = 1:length(Cluster_Vector(ii).SC)
    Silhouette_Total(end+1) = Cluster_Vector(1,ii).SC(a,1);
    end
end
end

%average the individual silhouette coefficients to get an overall
%silhouette score.
Silhouette_Score = mean(Silhouette_Total);

%Calculates silhouette score for each cluster.
for ii = 1:k
    Silhouette_Score_Cluster(1,ii) = mean(Cluster_Vector(ii).SC);
end

%% This sections sets up the data to create a silhouette coefficient graph

%Number of subjects in each cluster
for ii = 1:k
    Cluster_Subject_Total(ii) = length(Cluster_Subjects{1,ii});
end

%places silhouette coefficients for each cluster in a matrix to be used to
%graph
for ii = 1:k 
    Silhouette_Scores_Arranged{:,ii} = Cluster_Vector(ii).SC;
end

%silhouette coefficient graph
barh(1,Silhouette_Scores_Arranged{:,1},1.0,'r')
hold on
barh(2,Silhouette_Scores_Arranged{:,2},1.0,'b')
hold on
barh(3,Silhouette_Scores_Arranged{:,3},1.0,'g')
hold off
