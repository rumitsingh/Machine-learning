clear;
clc;
filename =('Lumbar_FlexExt_Kinematics.xlsx');%load file
Matrix = xlsread(filename);%read file
Matrix1 = xlsread(filename);

%The while loop and for loop below will calculate the dtw between the first
%vector and all the other vectors in the matrix. The while loop will then
%switch compare the distance between the 2nd vector and the other vectors
%in the matrix and this will be done for each vector in the matrix.
z = size(Matrix,2);
w = 1;
while z > 3 %While loop will stop when the size of the Matrix (z) will reach this threshold
n=1;%set n to start off at 1
%Calculates the dynamic time warping between each individual time-series
%compared to every other time-series in the matrix. This process is done
%for each time-series.
Hierarchical_DTW = zeros(size(Matrix,2),size(Matrix,2));
while n < size(Matrix,2)+1%While loop will run until all time-series are used
    for k=1:size(Matrix,2)%for loop will run until all time-series are used
    Hierarchical_DTW(n,k) = dtw(Matrix(:,n),Matrix(:,k)); 
    end
n = n+1;
end

%I use the sort function to rearrange the matrix from lowest to highest.
%This will set up the zeros into the first column
Hierarchical_DTW_Rearranged = zeros(size(Matrix,2),size(Matrix,2));
for k = 1:size(Matrix,2)
    Hierarchical_DTW_Rearranged(k,:) = sort(Hierarchical_DTW(k,:));
end

Hierarchical_DTW_Rearranged(:,1) = [];%This removes the first column containing all the zeros
Hierarchical_DTW_Min = Hierarchical_DTW_Rearranged(:,1);%Select the column that has the lowest DTW for each time-series.

Hierarchical_Cluster_DTW_Min{w} = Hierarchical_DTW_Min;%This variable will hold all the DTW min values for each run. 

%This will find the time-series that is closest to each time-series
Hierarchical_find = zeros(size(Matrix,2),1);
for k = 1:size(Matrix,2)
    Hierarchical_find(k,1) = find(Hierarchical_DTW(k,:)==Hierarchical_DTW_Min(k,1));
end

%This places all the time-series (numbers) that are closest to each time-series within the same cell 
Find_organized = cell(1,size(Matrix,2));
for k = 1:size(Matrix,2)
    Find_organized{k} = find(Hierarchical_find == k);
end 

%Some time-series had multiple other time series that are close to them,
%but we only want to pair each time-series up with one other time-series.
%Here, I find the DTW values for each time-series to figure out which
%time-series is the closest.
k=1;
Find_organize_DTW = cell(1,size(Matrix,2));
while k < size(Matrix,2)+1
    a = length(Find_organized{1,k});
        if isempty(Find_organized{k})
            Find_organize_DTW{k} = 0;   
        elseif a > 0
            for a = 1:length(Find_organized{1,k})
                b = Find_organized{1,k}(a,1);
                Find_organize_DTW{1,k}(a,1) = Hierarchical_DTW_Min(b,1);
            end
        end
    k=k+1;
end

%Finds the minimum DTW in each cell
Find_organize_DTW_min = cell(1,size(Matrix,2));
for k = 1:size(Matrix,2)
    Find_organize_DTW_min{k} = min(Find_organize_DTW{1,k});
end

%Finds which time-series is associated with the DTW value found previously
Find_organized_final = cell(1,size(Matrix,2));
for k = 1:size(Matrix,2)
    Find_organized_final{k} = find(Hierarchical_DTW_Min(:,1)==Find_organize_DTW_min{k});
end

%This deletes duplicate values in Find_organized_final variable
for k = 1:size(Matrix,2)
    if length(Find_organized_final{1,k}) ~= 1
        Find_organized_final{1,k} = Find_organized_final{1,k};
    elseif length(Find_organized_final{1,k}) == 1
        if isempty(Find_organized_final{1,Find_organized_final{1,k}})
            Find_organized_final{1,k} = Find_organized_final{1,k};
        elseif Find_organized{1,Find_organized_final{1,k}} > 0
            Find_organized_final{1,Find_organized_final{1,k}} = [];
        end
    end
end

a=1;
while a < size(Matrix,2)+1
    for k = 1:size(Matrix,2)
        b = length(Find_organized_final{1,a});
        c = length(Find_organized_final{1,k});
        if a == k
            Find_organized_final{1,a} = Find_organized_final{1,k};
        elseif a ~= k 
            if b == c
                if Find_organized_final{1,a} == Find_organized_final{1,k}
                Find_organized_final{1,k} = [];
                elseif Find_organized_final{1,a} ~= Find_organized_final{1,k}
                    Find_organized_final{1,k} = Find_organized_final{1,k};
                end
            end
        end
    end
    a=a+1;
end

%This adds values that are equal to k in cells that only have one other
%value. 
for k = 1:size(Matrix,2)
    if length(Find_organized_final{1,k}) == 1
        Find_organized_final{1,k}(2,1) = k;
    end
end

for k = 1:size(Matrix,2)
    if isempty(Find_organized_final{k})
        Find_organized_final{k} = Find_organized_final{k};
    elseif Find_organized_final{1,k}(1,1) > 0
        b = Find_organized_final{k}~=k;
        c = find(b==1);
        d = Find_organized_final{1,k}(c,1);
        Find_organized_final{1,d} = 0;
    end
end

%adds cell number in cells that are empty
%I do this because those values represent time-series that were not paired
%with another. 
for k = 1:size(Matrix,2)
    if isempty(Find_organized_final{1,k})
        Find_organized_final{1,k} = k;
    end
end

%delete cells that have a zero in them.
for k = 1:size(Matrix,2)
    if (Find_organized_final{1,k}) == 0
        Find_organized_final{1,k} = [];
    end
end

%Find all cells that are empty and move them to the end of the cell array
Hierarchical_Num = cell(1,size(Matrix,2));
for i = 1:size(Find_organized_final,1)
    j = ~cellfun('isempty',Find_organized_final(i,:));
    numNotEmpty = sum(j);
    Hierarchical_Num = repmat({missing},1,size(Find_organized_final,2));
    Hierarchical_Num(1,1:numNotEmpty) = Find_organized_final(i,j);
    Hierarchical_Num(i,:) = Hierarchical_Num;
end

%Delete all cells that are empty
Hierarchical_Num_Cut = cell(1,numNotEmpty);
for k = 1:numNotEmpty
    Hierarchical_Num_Cut{1,k} = Hierarchical_Num{1,k}; 
end

%This only runs on for the first time the while loop runs in order to
%transfer the values in Hierarchical_Num_Cut to Hierarchical_Num_Total
%This is so I can use Hierarchical_Num_Total for the rest of the runs
%through the while loop.
if w==1
    Hierarchical_Num_Total = Hierarchical_Num_Cut;
end

%This runs on all runs through the while loop after the first one.
%Hierarchical_Num_Total is a variable that holds all the subject numbers
%within each cluster. 
if w>1
    for b = 1:length(Hierarchical_Num_Cut)
        d = 0;
        for a = 1:length(Hierarchical_Num_Cut{1,b})
            for c = 1:length(Hierarchical_Num_Total{w-1,(Hierarchical_Num_Cut{1,b}(a,1))})
                if d == 0
                    Hierarchical_Num_Total{w,b}(1,1) = Hierarchical_Num_Total{w-1,(Hierarchical_Num_Cut{1,b}(a,1))}(1,1);
                elseif d > 0
                    Hierarchical_Num_Total{w,b}(end+1,1) = Hierarchical_Num_Total{w-1,(Hierarchical_Num_Cut{1,b}(a,1))}(c,1);
                end
                d = d+1;
            end
        end
    end
end
%Need to delete all contents of Hierarchical_Vector variable in order for
%it to run through all the while loops.
Hierarchical_Vectors = [];

%Time-series data from the Matrix are matched to their cluster.
for ii = 1:length(Hierarchical_Num_Total)
    for a = 1:length(Hierarchical_Num_Total{w,ii})
        Hierarchical_Vectors{1,ii}(:,a) = Matrix1(:,Hierarchical_Num_Total{w,ii}(a,1));
    end
end

%Need to delete all contents of Hierarchical_Vector_Sum variable in order for
%it to run through all the while loops.
Hierarchical_Vector_Sum = [];

%Need to delete all contents of Hierarchical_Vector_Ave variable in order for
%it to run through all the while loops.
Hierarchical_Vector_Ave = [];

%Calculates the sum of all the time-series within each cluster and averages
%them
for ii = 1:length(Hierarchical_Vectors)
    Hierarchical_Vector_Sum(:,ii) = zeros(length(Hierarchical_Vectors{1,1}(:,1)),1);
    for a = 1:size(Hierarchical_Vectors{1,ii},2)
        Hierarchical_Vector_Sum(:,ii) = Hierarchical_Vector_Sum(:,ii) + Hierarchical_Vectors{1,ii}(:,a);
    end
    Hierarchical_Vector_Ave(:,ii) = Hierarchical_Vector_Sum(:,ii)/size(Hierarchical_Vectors{1,ii},2);
end

%Holds all the subject number data for each while loop run.
Hierarchical_Num_Cluster{w} = Hierarchical_Num_Cut;
Matrix = Hierarchical_Vector_Ave;

%Finds the DTW min value that corresponds to each pairing of clusters.
for ii = 1:length(Hierarchical_Num_Cluster{w})
    Hierarchical_Cluster_DTW{1,w}{1,ii} = (Hierarchical_Cluster_DTW_Min{1,w}(Hierarchical_Num_Cluster{1,w}{1,ii}(1,1),1));
end

%Adds data from Hierarchical_Cluster_DTW to Hierarchical_Num_Cluster
Hierarchical_Num_Cluster{1,w}(2,:) = Hierarchical_Cluster_DTW{1,w}(1,:);

%This tracks the number of total clusters. This is used to terminate the
%while loop when the number of clusters reaches a certain level.
z = size(Matrix,2);

%This tells me how many times the while loop ran.
w = w+1;
end

for ii = 1:z
    Cluster_Subjects{1,ii} = Hierarchical_Num_Total{w-1,ii};
end

%figure
figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
hold on
for ii = 1:z
    plot(Matrix(:,ii),'LineWidth',3.0)
end
hold on
