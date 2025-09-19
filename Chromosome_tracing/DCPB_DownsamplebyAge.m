%% DCPB – 2 Jul 2024
% This script downsamples the chromosomes across all datasets based on age. 
% For each age group, it determines the minimum number of traces available 
% across all chromosome files and downsamples accordingly.

minRegions = 7;

load('/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/01_N2/N2_R15_30REG_AllAges_merge4_15_17.mat')

Chr_chosen = Chr(arrayfun(@(x) sum(x.r), Chr) >= minRegions);
Struct1 = Chr_chosen;

clearvars Chr
%mutant H3K9 

load('/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/02_GW_met2set25/GW_30REG_AllAges_merge10_14.mat')
Chr_chosen = Chr(arrayfun(@(x) sum(x.r), Chr) >= minRegions);
Struct2 = Chr_chosen;

clearvars Chr
%mutant condendsin 

load('/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/03_CondensinMT/MT_30REG_Allages_merge25_26_27.mat')
Chr_chosen = Chr(arrayfun(@(x) sum(x.r), Chr) >= minRegions);
Struct3 = Chr_chosen;
clearvars Chr

% Extract age fields
ages1 = [Struct1.Age];
ages2 = [Struct2.Age];
ages3 = [Struct3.Age];

uniqueAges1 = unique(ages1);
counts1= histc(ages1, uniqueAges1);

uniqueAges2 = unique(ages2);
counts2 = histc(ages2, uniqueAges2);

uniqueAges3 = unique(ages3);
counts3 = histc(ages3, uniqueAges3);

% To create a file that contains the ages and the minimum 
% number of counts among the three sets

% Combine all unique ages
allUniqueAges = unique([uniqueAges1, uniqueAges2, uniqueAges3]);

% Initialize arrays to store minimum counts
minCounts = zeros(size(allUniqueAges));

% Find the minimum count for each unique age across all sets
for i = 1:length(allUniqueAges)
    % Initialize with a large number
    minCount = Inf;
    
    % Check counts from each set and update minimum count
    if ismember(allUniqueAges(i), uniqueAges1)
        idx = find(uniqueAges1 == allUniqueAges(i));
        minCount = min(minCount, counts1(idx));
    end
    
    if ismember(allUniqueAges(i), uniqueAges2)
        idx = find(uniqueAges2 == allUniqueAges(i));
        minCount = min(minCount, counts2(idx));
    end
    
    if ismember(allUniqueAges(i), uniqueAges3)
        idx = find(uniqueAges3 == allUniqueAges(i));
        minCount = min(minCount, counts3(idx));
    end
    
    % Store the minimum count
    minCounts(i) = minCount;
end

%Remove the ages that are not in all data sets

%Initialize indices to remove
indicesToRemove = [];

% Iterate over allUniqueAges
for i = 1:length(allUniqueAges)
    % Initialize with a large number
    rev = Inf;

    % Check counts from each set and update minimum count
    if ~ismember(allUniqueAges(i), uniqueAges1)
        indicesToRemove(end+1) = i;
    elseif ~ismember(allUniqueAges(i), uniqueAges2)
        indicesToRemove(end+1) = i;
    elseif ~ismember(allUniqueAges(i), uniqueAges3)
        indicesToRemove(end+1) = i;
    end
end

% Remove elements from allUniqueAges and minCounts using logical indexing
allUniqueAges(indicesToRemove) = [];
minCounts(indicesToRemove) = [];


% Initialize an empty structure array to store results
ChrDS1 = struct([]);

% Iterate through each unique age

for i = 1:length(allUniqueAges)
    % Find indices in structArray where Age matches allUniqueAges(i)
    idx = find([Struct1.Age] == allUniqueAges(i));
    
    % Check if the count matches minCounts(i)
    if numel(idx) >= minCounts(i)
        % Randomly select minCounts(i) indices
        selectedIdx = idx(randperm(numel(idx), minCounts(i)));
        
        % Append selected entries to filteredStructArray
        ChrDS1 = [ChrDS1, Struct1(selectedIdx)]; 
    end
end

% Initialize an empty structure array to store results
ChrDS2 = struct([]);

% Iterate through each unique age
for i = 1:length(allUniqueAges)
    % Find indices in structArray where Age matches allUniqueAges(i)
    idx = find([Struct2.Age] == allUniqueAges(i));
    
    % Check if the count matches minCounts(i)
    if numel(idx) >= minCounts(i)
        % Randomly select minCounts(i) indices
        selectedIdx = idx(randperm(numel(idx), minCounts(i)));
        
        % Append selected entries to filteredStructArray
        ChrDS2 = [ChrDS2, Struct2(selectedIdx)]; 
    end
end

% Initialize an empty structure array to store results
ChrDS3 = struct([]);

% Iterate through each unique age
for i = 1:length(allUniqueAges)
    % Find indices in structArray where Age matches allUniqueAges(i)
    idx = find([Struct3.Age] == allUniqueAges(i));
    
    % Check if the count matches minCounts(i)
    if numel(idx) >= minCounts(i)
        % Randomly select minCounts(i) indices
        selectedIdx = idx(randperm(numel(idx), minCounts(i)));
        
        % Append selected entries to filteredStructArray
        ChrDS3 = [ChrDS3, Struct3(selectedIdx)]; 
    end
end

figure
hold on
histogram([ChrDS1.Age],10)
histogram([ChrDS2.Age],10)
histogram([ChrDS3.Age],10)
hold off

%% Downsamples using bin of ages 

clear all 
%close all

minRegions = 7;

% Load datasets
load('/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/01_N2/N2_R15_30REG_AllAges_merge4_15_17.mat')
Chr_chosen = Chr(arrayfun(@(x) sum(x.r), Chr) >= minRegions);
Struct1 = Chr_chosen;
clearvars Chr

load('/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/02_GW_met2set25/GW_30REG_AllAges_merge10_14.mat')
Chr_chosen = Chr(arrayfun(@(x) sum(x.r), Chr) >= minRegions);
Struct2 = Chr_chosen;
clearvars Chr

load('/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/03_CondensinMT/MT_30REG_Allages_merge25_26_27.mat')
Chr_chosen = Chr(arrayfun(@(x) sum(x.r), Chr) >= minRegions);
Struct3 = Chr_chosen;
clearvars Chr

% Define the age group bin size
binSize = 5;

% Bin ages by intervals of 4 years
binAges = @(age) floor(age/binSize)*binSize; 

% Extract and bin ages
ages1 = binAges([Struct1.Age]);
ages2 = binAges([Struct2.Age]);
ages3 = binAges([Struct3.Age]);

% Get unique binned ages and counts for each set
uniqueBins1 = unique(ages1);
counts1 = histc(ages1, uniqueBins1);

uniqueBins2 = unique(ages2);
counts2 = histc(ages2, uniqueBins2);

uniqueBins3 = unique(ages3);
counts3 = histc(ages3, uniqueBins3);

% Combine all unique bins
allUniqueBins = unique([uniqueBins1, uniqueBins2, uniqueBins3]);

% Initialize array for minimum counts across datasets for each age bin
minCounts = zeros(size(allUniqueBins));

% Find minimum counts for each bin
for i = 1:length(allUniqueBins)
    minCount = Inf;
    
    % Check counts from each dataset and update minimum count
    if ismember(allUniqueBins(i), uniqueBins1)
        idx = find(uniqueBins1 == allUniqueBins(i));
        minCount = min(minCount, counts1(idx));
    end
    
    if ismember(allUniqueBins(i), uniqueBins2)
        idx = find(uniqueBins2 == allUniqueBins(i));
        minCount = min(minCount, counts2(idx));
    end
    
    if ismember(allUniqueBins(i), uniqueBins3)
        idx = find(uniqueBins3 == allUniqueBins(i));
        minCount = min(minCount, counts3(idx));
    end
    
    % Store the minimum count for this bin
    minCounts(i) = minCount;
end

% Remove bins not present in all datasets
indicesToRemove = [];
for i = 1:length(allUniqueBins)
    if ~ismember(allUniqueBins(i), uniqueBins1) || ...
       ~ismember(allUniqueBins(i), uniqueBins2) || ...
       ~ismember(allUniqueBins(i), uniqueBins3)
        indicesToRemove(end+1) = i;
    end
end

allUniqueBins(indicesToRemove) = [];
minCounts(indicesToRemove) = [];

% Initialize structure arrays for downsampled data
ChrDS1 = struct([]);
ChrDS2 = struct([]);
ChrDS3 = struct([]);

% Downsample based on grouped age bins
for i = 1:length(allUniqueBins)
    % Downsample for Struct1
    idx1 = find(binAges([Struct1.Age]) == allUniqueBins(i));
    if numel(idx1) >= minCounts(i)
        selectedIdx1 = idx1(randperm(numel(idx1), minCounts(i)));
        ChrDS1 = [ChrDS1, Struct1(selectedIdx1)];
    end
    
    % Downsample for Struct2
    idx2 = find(binAges([Struct2.Age]) == allUniqueBins(i));
    if numel(idx2) >= minCounts(i)
        selectedIdx2 = idx2(randperm(numel(idx2), minCounts(i)));
        ChrDS2 = [ChrDS2, Struct2(selectedIdx2)];
    end
    
    % Downsample for Struct3
    idx3 = find(binAges([Struct3.Age]) == allUniqueBins(i));
    if numel(idx3) >= minCounts(i)
        selectedIdx3 = idx3(randperm(numel(idx3), minCounts(i)));
        ChrDS3 = [ChrDS3, Struct3(selectedIdx3)];
    end
end

% Plot histograms of downsampled age groups
figure
hold on
histogram([ChrDS1.Age], 'BinWidth', binSize, 'DisplayStyle', 'stairs')
histogram([ChrDS2.Age], 'BinWidth', binSize, 'DisplayStyle', 'stairs')
histogram([ChrDS3.Age], 'BinWidth', binSize, 'DisplayStyle', 'stairs')
hold off

% Plot histograms of downsampled age groups
figure
numBins = 15;
histogram([ChrDS1.Age] ,'FaceColor', [0 0.4470 0.7410], 'NumBins', numBins);
hold on;
histogram([ChrDS2.Age], 'FaceColor', [0.6350 0.0780 0.1840], 'FaceAlpha', 0.3, 'NumBins', numBins);
histogram([ChrDS3.Age], 'FaceColor', [0.1660 0.5740 0.1880], 'FaceAlpha', 0.3, 'NumBins', numBins);
legend('Wild Type', '(-)H3K9me', '(-)Condensin I');
xlabel('Age');
ylabel('Number of Traces');
hold off;


%% Generate matrices

RegionStart = 1 

% load ax values for axis labels
load('1to40_ce10_axValues.mat')
axValues  = axValues(14:40);

SkipTADs = sort([28, 29, 40], 'descend') -13;
varList = {'ChrDS1', 'ChrDS2', 'ChrDS3'};

for l= 1:length(varList)
    TotalNumTADs= 30;
    varName = varList{l};
    Chr_chosen = eval(varName);

    Mean = zeros(TotalNumTADs);

    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            DisList = zeros(1, length(Chr_chosen)) - 1;
            for k = 1:length(Chr_chosen)

                if Chr_chosen(k).r(i) == 1 && Chr_chosen(k).r(j) == 1 %%&& Chr_chosen(k).Age >= AgeStart && Chr_chosen(k).Age <= AgeStop %silence the age range if you DO NOT need it
                    if ((Chr_chosen(k).x(i)-Chr_chosen(k).x(j))^2+(Chr_chosen(k).y(i)-Chr_chosen(k).y(j))^2+(Chr_chosen(k).z(i)-Chr_chosen(k).z(j))^2)^0.5 <=6
                        DisList = [DisList ((Chr_chosen(k).x(i)-Chr_chosen(k).x(j))^2+(Chr_chosen(k).y(i)-Chr_chosen(k).y(j))^2+(Chr_chosen(k).z(i)-Chr_chosen(k).z(j))^2)^0.5];
                    end
                end

                % Confirm distance ≤ 6 and ensures we don't have NaN's
                % if DisList <= 6 && ~isnan(DisList)
                %     DisList(k) = DisList;
                % end
            end
            DisList(DisList == -1) = [];

            Mean(i, j) = median(DisList); %ATTENTION:i changed for the median
            DisListAll{i,j} = DisList;
        end
    end

    for i = 1:length(SkipTADs)
        Mean(SkipTADs(i), :) = [];
        Mean(:, SkipTADs(i)) = [];
    end

    TotalNumTADs = length(Mean);

    figure
    imagesc(Mean);
    set(gca, 'colormap', hot)
    RedBlue
    caxis([0, 1.5])
    % Uncomment the following line to preserve the X-limits of the axes
    xlim([RegionStart - 0.5 TotalNumTADs + 0.5]); %28.5
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim([RegionStart - 0.5 TotalNumTADs + 0.5]);
    yticks(RegionStart:5:TotalNumTADs)
    yticklabels({axValues(RegionStart:5:TotalNumTADs)})
    xticks(RegionStart:5:TotalNumTADs)
    xticklabels({axValues(RegionStart:5:TotalNumTADs)})
    xtickangle(90)
    colorbar;
    axis square
    title(['Median of Distance Matrices - n=' num2str(length(Chr_chosen))])
    %set(gcf,'units','centimeters','position',[10,10,30,30])
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'FontName', 'Arial');

end
