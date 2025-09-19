%% REG30 Clusters Trace Boundary Calculation Clusters - 2024 
% This script was written by D.C.P.B. and incorporates improvements from T. Verheijen.

% This script computes the probabilities of the locations of boundaries
% within single cells. The insulation is calculated like in crane et al.,
% 2015 with the squares along the diagonal. The prevalence of the boundary
% is calculated similarly to Chen et al., 2021 when the data is filled
% linearly and after the insulation (like Crane), the number of boundaries
% is divided in the among of total single traces.


close all
clearvars -except y_fit

% DECLARE if you would like to save the figures by setting savefigs = 1,
% else 0
savefigs = 1;
addpath('/scicore/home/mangos/pulido0000/Tracing/BasicFunctions')

% load ax values for axis labels
load('1to40_ce10_axValues.mat')
axValues  = axValues(1:13);

%load the average insualtion just for the plot 
%wt
%load('/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/01_N2/F_insulation_median_WT_2to140__min_7_TH_0.2_TadNum_27/InsulationProfileNorm_2_Delta_1_Profile_TH_0.2_Reg_1to30.mat')
%met-2
%load('/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/02_GW_met2set25/insulation_median_GW_2to140__min_7_TH_0.1_TadNum_27/InsulationProfileNorm_2_Delta_1_Profile_TH_0.1_Reg_1to30.mat')
%condensin 
%load('/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/03_CondensinMT/insulation_median_MT_2to140__min_7_TH_0.1_TadNum_27/InsulationProfileNorm_2_Delta_1_Profile_TH_0.1_Reg_1to30.mat')

load('/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/CORG/N2_insulation_median_R15_2to140__min_1_TH_0.1_TadNum_13/InsulationProfileNorm_2_Delta_1_Profile_TH_0.1_Reg_1to13.mat')

InsulationAverage = InsProfNorm;

%Input this if you want to normalize the insulation to the average
%MeanInsulationAverage=  0.9515; %WT
%MeanInsulationAverage=  0.9486; %MT
%MeanInsulationAverage=  0.9790; %MT


% load the Y_Fit of the mean if you are using it
% load('old_data/Min14_2to40cells_allchr_Y_Fit.mat');

% name of the cluster folder from the current path
foldername = 'Res 0.8';
Resolution = extractAfter(foldername, strlength(foldername) - 3);
files_in_dir = dir([foldername '/*.mat']);
num_files = length(files_in_dir);

% Create colormap for Boundary Strength
StrengthColor = flipud(parula(450));
StrengthColor = StrengthColor(80:430, :);

fullParula = flipud(turbo(150));
% Find the approximate starting index for yellow and end index for dark blue
yellowStartIndex = 2;   % Approximate start index for yellow in parula
darkBlueEndIndex = 150;  % End index for the darkest blue
% Extract only the yellow to dark blue range
yellowToBlueRange = fullParula(yellowStartIndex:darkBlueEndIndex, :);
% Interpolate to get 100 colors with strong transitions
StrengthColor = interp1(linspace(0, 1, size(yellowToBlueRange, 1)), yellowToBlueRange, linspace(0, 1, 100));


% Define constants:
RegionStart = 1;
minRegions = 5;
Boundary_TH = 0.2;
TotalNumTADs = 13;
SkipTADs = sort([], 'descend');
SkipTADs = SkipTADs(SkipTADs < TotalNumTADs);
numRegions = TotalNumTADs - RegionStart - length(SkipTADs) + 1;

name = 'Cluster_Insulation_Res05_CORG';
dirname = ([num2str(name) '_min_' num2str(minRegions) '_TH_' ...
            num2str(Boundary_TH) '_Reg_' num2str(RegionStart) 'to' ...
            num2str(TotalNumTADs)]);
if ~isfolder(dirname) && savefigs == 1
    mkdir(dirname);
end

% Select quality Cells with high region counts
% Chr_chosen = Chr(arrayfun(@(x) sum(x.r), Chr) >= minRegions);
% Chr = Chr_chosen;
% generate_region_histogram(Chr_chosen);
% saveas(gcf, [dirname '/1to14_Num_Regions_Distribution_singles.fig'])

% Used to normalize matrices with y_fit of the average

fid = fopen('30reg_Probes_startends2829_40_Ce10.txt', 'r');
tline = fgetl(fid);
linenum = 1;
DomainStart = zeros(100, 1) - 1;
DomainEnd = zeros(100, 1) - 1;
DomainCenter = zeros(100, 1) - 1;
while ischar(tline)
    C = strsplit(tline);
    DomainStart(linenum) = str2double(C{2});
    DomainEnd(linenum) = str2double(C{3});
    DomainCenter(linenum) = (DomainStart(linenum) + DomainEnd(linenum))/2.0;
    tline = fgetl(fid);
    linenum = linenum + 1;
end
fclose(fid);

DomainStart = DomainStart(1:numRegions);
DomainEnd = DomainEnd(1:numRegions);
DomainCenter = DomainCenter(1:numRegions);
DomainCenter = DomainCenter./1000000; %convert from bp to Mb

DomainSize = length(DomainCenter);
GenomicDis = zeros(DomainSize);

for i = 1:DomainSize
    for j = 1:DomainSize
        GenomicDis(i, j) = abs(DomainCenter(i) - DomainCenter(j));
    end
end
GenomicDis_pxl = GenomicDis(tril(true(size(GenomicDis)), -1));

InsBoxSize = 2;
DeltaBoxSize = 1;
InsProfRange = InsBoxSize + 1:numRegions - InsBoxSize;
DeltaRange = InsBoxSize + DeltaBoxSize + 1:numRegions - ...
    InsBoxSize - DeltaBoxSize;

AllStrength = cell(num_files, 1);
AllBoundaries = cell(num_files, 1);

% tracesmissing = zeros(length(Chr_chosen), 1) - 1;

Cluster_Means = cell(num_files, 1);
InsProfNorm = zeros(num_files, length(InsProfRange));

numskipped = 0;
while sum(RegionStart == SkipTADs) ~= 0
    SkipTADs(end) = [];
    RegionStart = RegionStart + 1;
    numskipped = numskipped + 1;
end

%% Main Loop: Compute Boundary positions and strength
for k = 1:num_files
    numRegions = TotalNumTADs - RegionStart + 1;
    load([foldername '/' files_in_dir(k).name])
    
    %ONLY THIS FOR THE CoClustering
    % % Load the file into a temporary structure for the split ones 
    % temp = load(fullfile(foldername, files_in_dir(k).name));
    % % Extract the variable dynamically and assign it to 'ChrClus'
    % varName = fieldnames(temp); % Get the name of the loaded variable
    % ChrClus = temp.(varName{1}); % Assign the loaded variable to 'ChrClus'
    % clear temp;

    Cluster_Means{k} = computeOverallMean(ChrClus, SkipTADs - RegionStart + 1, numRegions);
    
    %just for the met-2, silence for the other strains, remove region 11
    %for the analysis just becasue low hyb numbers      
   %Cluster_Means{k}(11, :) = NaN;
   %Cluster_Means{k}(:, 11) = NaN;
    % 
    Mean_plot = Cluster_Means{k};
    numRegions = length(Cluster_Means{k});

    % Fill the Mean matrix with Linear Interpolation
    Mean_filtered = fillmissing(Cluster_Means{k}, 'linear'); % Fill Rows
    Mean_filtered = fillmissing(Mean_filtered, 'linear', 2); % Fill Columns
    Mean_filtered(11,11) = 0;

    % negatives = find(Mean_filtered < 0);
    % if ~isempty(negatives)
    %     Mean_filtered(negatives) = MeanofMeans(negatives);
    % end

    % If the Genomic Normalization is being calculated during each loop,
    % uncomment this section and comment out the Mean_adjust loop below
    % 
    % Dis_pxl = Mean_filtered(tril(true(size(Mean_filtered)), -1));
    % R = corrcoef(GenomicDis_pxl, Dis_pxl);
    % x = log(GenomicDis_pxl);
    % y = log(Dis_pxl);
    % X = [ones(size(x)), x];
    % %     [b, bint] = regress(y, X);
    % [b, bint] = regress(y, X);
    % x = GenomicDis_pxl';
    % y_fit = exp(b(1)) * x.^b(2);
    % norms = Dis_pxl./y_fit';
    % Mean_adjust = squareform(norms) + eye(numRegions);

    % Perform Genomic Normalization based with read in y_fit. Comment this
    % loop out if doing individual Genomic Normalization
    % 
    % Mean_adjust = Mean_filtered;
    % n = 1;
    % for i = 1:numRegions
    %     for j = 1:numRegions
    %         if i~=j
    %             Mean_adjust(i,j) = Mean_filtered(i,j)/y_fit(n);
    %             n = n + 1;
    %         else
    %             Mean_adjust(i,j) = 1;
    %         end
    %     end
    % end

    % Perform Genomic Normalization based with read in y_fit. Comment this
    % loop out if doing individual Genomic Normalization

    Mean_adjust = Mean_filtered;
    Dis_pxl = Mean_filtered(tril(true(size(Mean_filtered)), -1));
    norms = Dis_pxl./y_fit';
    Mean_adjust = squareform(norms) + eye(numRegions);
    Cluster_Means_adjust{k,1} = Mean_adjust; 

    % Find the mean value in an InsBoxSize * InsBoxSize matrix off of the
    % principal diagonal of the Mean_adjust matrix
    InsProfile = zeros(numRegions - 2 * InsBoxSize, 1);
    % We include -2 * InsBoxSize because we lose an InsBoxSize amount of
    % possible means off of either side of the matrix
    for i = 1:length(InsProfile)
        l = i + InsBoxSize;
        InsProfile(i) = 1.0/mean(Mean_adjust(i:l-1, ...
            l + 1:l + InsBoxSize), 'all');
    end

    % Normalize the Insulation Profiles
    InsProfNorm(k, :) = log2(InsProfile/mean(InsProfile));

    %Normalize the insulation to the average (Silence the one before)
    %InsProfNorm(k, :) = log2(InsProfile/MeanInsulationAverage);

    % Find the change in insulation values
    DeltaIns = zeros(length(InsProfNorm(k, :)) - 2 * DeltaBoxSize, 1);
    % We include -2 * DeltaBoxSize, because we cannot use the two
    % DeltaBoxSize means on either end of the InsProfNorm array
    for i = 1:length(DeltaIns)
        DeltaIns(i) = mean(InsProfNorm(k, i:i + DeltaBoxSize - 1)) - ...
            mean(InsProfNorm(k, i + DeltaBoxSize + 1:i + 2 * DeltaBoxSize));
    end

    % Find peaks and valleys of the insulation
    [InsMax, InsMaxIdx] = findpeaks(InsProfNorm(k, :));
    [InsMin, InsMinIdx] = findpeaks(-InsProfNorm(k, :));
    InsMin = -InsMin;
    InsMaxIdx = InsMaxIdx + InsBoxSize;
    InsMinIdx = InsMinIdx + InsBoxSize;

    InsPeaks = struct('InsMaxIdx', num2cell(InsMaxIdx), ...
        'InsMax', num2cell(InsMax));
    InsValleys = struct('InsMinIdx', num2cell(InsMinIdx), ...
        'InsMin', num2cell(InsMin));

    % Find peaks and valleys of the detla
    [DeltaMax, DeltaMaxIdx] = findpeaks(DeltaIns);
    [DeltaMin, DeltaMinIdx] = findpeaks(-DeltaIns);
    DeltaMin = -DeltaMin;
    DeltaMaxIdx = DeltaMaxIdx + InsBoxSize + DeltaBoxSize;
    DeltaMinIdx = DeltaMinIdx + InsBoxSize + DeltaBoxSize;

    DeltaPeaks = struct('DeltaMaxIdx', num2cell(DeltaMaxIdx), ...
        'DeltaMax', num2cell(DeltaMax));
    DeltaValleys= struct('DeltaMinIdx', num2cell(DeltaMinIdx), ...
        'DeltaMin', num2cell(DeltaMin));

    % This line of code is doing something interesting. First, we convert
    % yData to a logical array of 1s and 0s where a 1 indicates that it is
    % positive, and a 0 indicates that it is negative. Then we filter
    % through the results to find these crossover points.
    BoundaryPositions = sort(strfind(DeltaIns.' >= 0, [1 0]));

    % Here we check to find wich side of the zerocrossover is closer to
    % zero.
    for i = 1:length(BoundaryPositions)
        if BoundaryPositions(i) == length(DeltaRange)
            % Ensure that we won't extand past array bounds
            continue;
        end
        if abs(DeltaIns(BoundaryPositions(i))) > abs(DeltaIns(BoundaryPositions(i) + 1))
            BoundaryPositions(i) = BoundaryPositions(i) + 1;
        end
    end

    % Find the closest Delta Peaks and Valleys for boundary strengths
    % Find the closest Insulation Valleys for crossover point comparison
    ClosestDeltaPeakIdx = zeros(length(BoundaryPositions), 1);
    ClosestDeltaValleyIdx = zeros(length(BoundaryPositions), 1);
    ClosestInsValley = zeros(length(BoundaryPositions), 1);
    for i = 1:length(BoundaryPositions)
        % Find largest peak index value smaller than the boundary position
        if max(DeltaMaxIdx(DeltaMaxIdx <= DeltaRange(BoundaryPositions(...
                i)))) >= DeltaRange(BoundaryPositions(i)) - 3
            % Include -3 to remove small peaks that hide the next true peak
            % in the isulation

            ClosestDeltaPeakIdx(i) = find(DeltaMaxIdx == max(DeltaMaxIdx(...
                DeltaMaxIdx <= DeltaRange(BoundaryPositions(i)))));
        else
            % Flag as non-accepted peak
            ClosestDeltaPeakIdx(i) = -1;
        end

        % Find smallest peak index value greater than the boundary position
        if min(DeltaMinIdx(DeltaMinIdx >= DeltaRange(BoundaryPositions(...
                i)))) <= DeltaRange(BoundaryPositions(i)) + 4
            % Include +4 to remove small valleys that hide the next true
            % valley in the insulation
            ClosestDeltaValleyIdx(i) = find(DeltaMinIdx == min(DeltaMinIdx(...
                DeltaMinIdx >= DeltaRange(BoundaryPositions(i)))));
        else
            % Flag as non-accepted valley
            ClosestDeltaValleyIdx(i) = -1;
        end

        % Find indices where there is an insulation valley within ± 1 of a
        % zero crossover
        CloseInsVal = InsMinIdx(abs(InsMinIdx - ...
            DeltaRange(BoundaryPositions(i))) < 2);
        if ~isempty(CloseInsVal)
            % Default to the smallest such valley
            ClosestInsValley(i) = min(CloseInsVal);
        else
            % Flag as non-accepted position
            ClosestInsValley(i) = -1;
        end
    end

    % Compute boundary strength or reject invalid crossovers
    Boundaries = zeros(length(BoundaryPositions), 2);
    for i = 1:length(BoundaryPositions)
        if ClosestInsValley(i) == -1
            % Here the Zero crossover point is more than 2 positions
            % away from the Insulation valley, so we can skip this
            % boundary position
            Boundaries(i, :) = -1;
            continue;
        end
        % if either the nearest Delta Peak or Valley is -1, we know one is
        % invalid and thus we continue
        if ClosestDeltaPeakIdx(i) == -1 || ClosestDeltaValleyIdx(i) == -1
            Boundaries(i, :) = -1;
            continue;
        else
            % Both the nearest valley and peak are valid, so find strength
            strength = DeltaMax(ClosestDeltaPeakIdx(i)) - ...
                DeltaMin(ClosestDeltaValleyIdx(i));
            Boundaries(i, 1) = DeltaRange(BoundaryPositions(i));
            Boundaries(i, 2) = strength;
        end
    end

    % Remove invalid boundaries
    Boundaries = Boundaries(Boundaries(:, 1) ~= -1, :);

    % If no valid boundaries, skip to next trace
    %         if isempty(Boundaries)
    %             continue;
    %         end

    % Filter for boundaries above the boundary threshold
    Boundaries = Boundaries(Boundaries(:, 2) > Boundary_TH, :);

    [~, ind, ~] = unique(Boundaries(:, 2));
    Boundaries = Boundaries(sort(ind, 'ascend'), :); % Ensure lo-high order

    %         if k <= length(Chr_chosen)
    AllBoundaries{k} =  Boundaries(:, 1);
    AllStrength{k} =  Boundaries(:, 2);
    %         end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Draw figures for this trace

    t = draw_figs(Cluster_Means{k}, axValues, k, Mean_plot, ...
        RegionStart, TotalNumTADs, numRegions, Boundary_TH, ...
        InsProfRange, InsProfNorm(k, :), ...
        DeltaRange, DeltaIns,InsulationAverage,...% DomainCenter,
        Boundaries, StrengthColor, dirname, InsBoxSize, DeltaBoxSize, ...
        savefigs, length(SkipTADs), numskipped, ChrClus);
     % 
     % u = draw_figsNorm(Cluster_Means_adjust{k}, axValues, k, Mean_adjust, ...
     %    RegionStart, TotalNumTADs, numRegions, Boundary_TH, ...
     %    InsProfRange, InsProfNorm(k, :), ...
     %    DeltaRange, DeltaIns,...% DomainCenter,
     %    Boundaries, StrengthColor, dirname, InsBoxSize, DeltaBoxSize, ...
     %    savefigs, length(SkipTADs), numskipped, ChrClus);

    %close all
end % End for 1:length(Chr_chosen)

AllBoundaries = AllBoundaries(~cellfun('isempty', AllBoundaries));
AllStrength = AllStrength(~cellfun('isempty', AllStrength));
% save([dirname '/1to14_AllBoundaries_TH_' num2str(Boundary_TH) '_min_' num2str(minRegions) '.mat'], 'AllBoundaries');
% save([dirname '/1to14_AllStrength_TH_' num2str(Boundary_TH) '_min_' num2str(minRegions) '.mat'], 'AllStrength');

total_occurrences = zeros(numRegions, 1);
for i = 1:length(total_occurrences)
    Counts = cellfun(@(x) find(x == i), AllBoundaries, 'UniformOutput', 0);
    total_occurrences(i) = numel([Counts{:}]);
end

% The new line is uncommented above the old line
total_occurrences = total_occurrences/length(AllBoundaries);


%% Boundary Probability plot

if ~isempty(AllBoundaries)
    figure
    hold on
    ydata = 1:numRegions; %this is the x data
    sz = numRegions;
    r= scatter(ydata, total_occurrences(:), sz, 'filled');
    rr= plot(total_occurrences(:), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
    yline(mean(nonzeros(total_occurrences)), '--', 'Mean', 'LineWidth', 3);
    grid on
    yline(mean(nonzeros(total_occurrences)) + std(nonzeros(total_occurrences)), ':', 'LineWidth', 2);
    yline(mean(nonzeros(total_occurrences)) - std(nonzeros(total_occurrences)), ':', 'LineWidth', 2);
    xlim([1 numRegions]);
    xticks(1:numRegions);
    xticklabels({axValues(RegionStart - numskipped:TotalNumTADs - length(SkipTADs) - numskipped)});
    %xticklabels({round(DomainCenter(1,RegionStart:TotalNumTADs),1)})
    xtickangle(90);
    ylim([0 1])
    title('Boundary probability from single cells ', 'FontSize', 12, 'fontname', 'Arial', 'VerticalAlignment', 'bottom');
    ylabel('Boundary probability', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold');
    xlabel('Region', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold')%xL = xlim;
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'FontName', 'Arial');
    hold off
    
    if savefigs == 1
        saveas(gcf, [dirname '/Boundary_Probability_Insulation_' ...
                    num2str(InsBoxSize) '_Delta_' num2str(DeltaBoxSize) ...
                    '_Cluster_Res_' Resolution '.fig'])
    end
end

%close all

%% Create All Insulations figure

fig = figure;
hold on
ydata = 1:numRegions;
% plot(ydata, InsProfNorm{:}, 'LineWidth', 1)
% plot(ydata, mean(InsProfNorm), ':','LineWidth', 1.5, 'color', [0 0 0])
plot(InsProfRange, mean(InsProfNorm), 'LineWidth', 1.5, 'color', ...
                                    [0 0 0], 'HandleVisibility', 'on', 'DisplayName', 'Mean')

cmap = turbo(length(files_in_dir) + 4);
legnames = strings(length(files_in_dir) + 1, 1);
legnames(1) = 'Means';
for k = 1:length(files_in_dir)
    if k <= length(files_in_dir)/2
        plot(InsProfRange, InsProfNorm(k,: ), 'LineWidth', 1.5, 'color',...
                                    cmap(k + 1, :), 'HandleVisibility', 'on', 'DisplayName', ['Clus' num2str(k - 1)])
    else
        plot(InsProfRange, InsProfNorm(k,: ), 'LineWidth', 1.5, 'color',...
                                cmap(3 + k, :), 'HandleVisibility', 'on', 'DisplayName', ['Clus' num2str(k - 1)])
    end
    legnames(k + 1) = ['Cluster ' num2str(k - 1)];
end
yticks(-1.1:0.2:1.1)
ylim([-1.1 1.1])
xlim([1 numRegions])
xticks(1:numRegions)
xticklabels({axValues(RegionStart - numskipped:TotalNumTADs - length(SkipTADs) - numskipped)})
xtickangle(90)
title(['Cluster Insulation Profiles'])
ylabel('Insulation Score')
xlabel('Region')
grid on
plot(InsProfRange, mean(InsProfNorm), 'LineWidth', 8, 'color', ...
    [0.7 0.7 0.7], 'HandleVisibility', 'on', 'DisplayName', 'Mean_shadow')
ax = gca;
set(gca, 'FontName', 'Arial', 'FontSize', 12)
hold off
[ll, lgd, gg, ~] = legend(legnames);
hold on
lgd_lines = findobj(lgd, 'type', 'line');
lgd_lbls = findobj(lgd, 'type', 'text');
set(lgd_lbls, 'FontName', 'Arial', 'FontSize', 24)
set(lgd_lines, 'LineWidth', 2);
hold off
ax.Children = ax.Children(end:-1:1);

if savefigs == 1
    saveas(fig, [dirname '/All_Insulation_Profiles_Insulation_' ...
                num2str(InsBoxSize) '_Res_' Resolution '_Clusters.fig']);
end

%close

%% Create SEM figures

figure
stdshade(InsProfNorm(:, :), 0.1, [0 102 180]/255); %sem shade in the fcn code
% xlim([RegionStart TotalNumTADs] - 2)
% xticks(RegionStart - 2:TotalNumTADs - 2)
% xticklabels({round(DomainCenter(RegionStart:end), 1)})
xlim([1 numRegions] - InsBoxSize)
xticks(1 - InsBoxSize:numRegions - InsBoxSize)
xticklabels({axValues(RegionStart - numskipped:TotalNumTADs - length(SkipTADs) - numskipped)})
xtickangle(90)
yticks(-1.1:0.2:1.1)
ylim([-1.1 1.1])
title('Cluster Insulation SEM')
grid on
ylabel('Insulation Score')
xlabel('Region')
set(gca, 'FontName', 'Arial', 'FontSize', 15)

if savefigs == 1
    saveas(gcf, [dirname '/All_Insulation_Profiles_SEMshade_Res_' ...
                                            Resolution '_Clusters.fig']);
end


figure
stdshade(InsProfNorm(:, :), 0.1, [0 102 180]/255); %sem shade in the fcn code
hold on
boxplot(InsProfNorm(:, :))
% xlim([RegionStart TotalNumTADs] - 2)
% xticks(RegionStart - 2:TotalNumTADs - 2)
% xticklabels({round(DomainCenter(RegionStart:end), 1)})
xlim([1 numRegions] - InsBoxSize)
xticks(1 - InsBoxSize:numRegions - InsBoxSize)
xticklabels({axValues(RegionStart - numskipped:TotalNumTADs - length(SkipTADs) - numskipped)})
xtickangle(90)
yticks(-1.1:0.2:1.1)
ylim([-1.1 1.1])
title('Cluster Insulation SEM with Regional Boxplots')
grid on
%ylim([-0.5 0.5])
ylabel('Insulation Score')
xlabel('Region')
set(gca, 'FontName', 'Arial', 'FontSize', 35)

if savefigs == 1
    saveas(gcf, [dirname '/All_Insulation_Profiles_SEMshade_Res_' ... 
                                    Resolution '_Clusters_xBoxPlot.fig']);
end

%close all
StrengthColor

%% Declare Functions

function generate_region_histogram(a)
    figure
    histogram(arrayfun(@(x) sum(x.r), a))
    xlabel('Number of regions detected')
    ylabel('Number of traces')
    ColorMap = load('~/BasicFunctions/RedBlue.txt');
    colormap(ColorMap/255);
end

function Mean = computeOverallMean(C, SkipTADs, numRegions)
    Mean = zeros(numRegions);
    
    for i = 1:numRegions
        for j = 1:numRegions
            DisList = zeros(1, length(C)) - 1;
            for k = 1:length(C)
%                 disp([k, i, j])
                if ~C(k).r(i) || ~C(k).r(j)
                    continue;
                end
    
                dist = distance(C(k), i, j);
                % Confirm distance ≤ 6 and ensures we don't have NaN's
                if dist <= 6 && ~isnan(dist)
                    DisList(k) = dist;
                end
            end
            DisList(DisList == -1) = [];
    
            Mean(i, j) = median(DisList); %%NOTE THIS IS CHANGE FOR THE MEADIAN
        end
    end
    
    % SkipTADs = SkipTADs(SkipTADs > RegionStart) - RegionStart + 1;
    Mean = remove_TADs(Mean, SkipTADs);
end

function Dist = computeDistMat(C, SkipTADs, NumRegions)
    Dist = zeros(NumRegions);

    for i = 1:NumRegions
        for j = 1:NumRegions
            if C.r(i) == 1 && C.r(j) == 1
                Dist(i,j) = distance(C, i, j);
            else
                Dist(i,j) = NaN;
            end
        end
    end

    Dist = remove_TADs(Dist, SkipTADs);
end

function mat = remove_TADs(mat, skips)
    for i = 1:length(skips)
        mat(skips(i), :) = [];
        mat(:, skips(i)) = [];
    end
end

function d = distance(vect, i, j)
% Compute distance given a 3-dimensional vector as input
xsq = (vect.x(i) - vect.x(j))^2;
ysq = (vect.y(i) - vect.y(j))^2;
zsq = (vect.z(i) - vect.z(j))^2;
d = sqrt(xsq + ysq + zsq);
end

function t = draw_figs(Mean, axValues, k, Mean_plot, RegionStart, ...
                                TotalNumTADs,numRegions, Boundary_TH, ...
                                InsProfRange, InsulationProfNorm, ...
                                DeltaRange, DeltaInsulation, InsulationAverage,... % DomainCenter, ...
                                Boundaries, StrengthColor, dirname, ...
                                InsBoxSize, DeltaBoxSize, s, numskip, numskipped,ChrClus)
    
    figure
    t = tiledlayout(5, 1, 'TileSpacing', 'compact');
    
    nexttile(1, [3 1])
    RedBlue
    h= imagesc(Mean_plot);
    
    h.AlphaData = ~isnan(Mean_plot); %to make the isnan values another color 
    set(gca, 'Color', [0.6 0.6 0.8])
    
    % Uncomment the following line to preserve the X-limits of the axes
    xlim([0.5 numRegions + 0.5]); %28.5
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim([0.5 numRegions + 0.5]);
    yticks(1:2:numRegions + 1);
    yticklabels({axValues(RegionStart - numskipped:2:TotalNumTADs - numskip - numskipped)})
    xticks(1:2:numRegions + 1);
    xticklabels({axValues(RegionStart - numskipped:2:TotalNumTADs - numskip - numskipped)})
    xtickangle(90);
    colorbar;
    axis square
    caxis([0, 1.5]);
    title(['Median distance Cluster ' num2str(k - 1) ' n=' num2str(length(ChrClus))]);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'FontName', 'Arial');

    nexttile(4, [2 1])
    %nexttile(10,[2 1])
    hold on
    %     ydata = RegionStart:TotalNumTADs;
    %h = yline(0, 'k-', 'LineWidth', 1);
    %         area(InsulationProfilesNorm(p+1, :), 'LineWidth', 1, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7]);
    %         r = plot(InsulationProfilesNorm(p+1, :), 'LineWidth', 1, 'Color', [0 0 0]);
    %         rr = plot(DeltaInsulation(p+1, :), 'LineWidth', 1, 'Color', [1 0 0]);
    area(InsProfRange, InsulationProfNorm(:), 'LineWidth', 1, ...
        'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7]);
    Ins = plot(InsProfRange, InsulationProfNorm(:), 'LineWidth', ...
        1, 'Color', [0 0 0]);

    rDelta = plot(DeltaRange, DeltaInsulation(:), 'LineWidth', 1, ...
            'Color', [1 0 0]);r = plot(DeltaRange, DeltaInsulation(:), 'LineWidth', 1, ...
         'Color', [1 0 0]);
    plot(InsProfRange, InsulationAverage, 'LineWidth', 1.5, 'color', ...
                                    [0 0.4470 0.7410], 'HandleVisibility', 'on', 'DisplayName', 'Mean')
    xlim([0.5 numRegions + 0.5]);
    xticks(1:2:numRegions);
    xticklabels({axValues(RegionStart - numskipped:2:TotalNumTADs - numskip - numskipped)});
    % xticklabels({round(DomainCenter(1,RegionStart:TotalNumTADs),1)});
    xtickangle(90)
    title('Insulation(Black) and delta(Red) profile');
    %grid on
    ylabel('score');
    xlabel('region');
    ylim([-1.1 1.1]);
    % yL = ylim;
    if ~isempty(Boundaries) %if there are not boundaries this won't run
        for i = 1:length(Boundaries(:, 1))
            if Boundaries(i, 2) > 1.0
                %line([Boundaries(i,1) Boundaries(i,1)],[0.52 0.55], yL, 'color',[0.6 0.6 0.6], 'linewidth', 3 );
%                 rrr = line([Boundaries(i, 1) Boundaries(i, 1)], [1.00 1.03], yL, 'color', StrengthColor(round(Boundaries(i, 2) * 100), :, :), 'linewidth', 3);
                bound = scatter(Boundaries(i, 1), 1.0, 50, StrengthColor(round(Boundaries(i, 2) * 10), :), 'v', 'filled');
            elseif Boundaries(i, 2) > Boundary_TH % this is 0.1 th of crane.
%                 rrr = line([Boundaries(i, 1) Boundaries(i, 1)], [Boundaries(i, 2) (Boundaries(i, 2) + 0.03)], yL, 'color', StrengthColor(round(Boundaries(i, 2) * 100), :), 'linewidth', 6);
                bound = scatter(Boundaries(i, 1), Boundaries(i, 2), [], StrengthColor(round(Boundaries(i, 2) * 100), :), 'v', 'filled');
            end%x-axis
        end%x-axis
        legend([Ins, bound], 'Insulation', 'Boundary', 'FontSize', 12, 'fontname', 'Arial', 'Location', 'southeast');
    else
        legend(Ins, 'Insulation', 'FontSize', 12, 'fontname', 'Arial', 'Location', 'southeast');
    end
    
    title('Insulation and delta profile', 'FontSize', 12, 'fontname', 'Arial', 'VerticalAlignment', 'bottom');
%     legend([r, rr], 'Insulation', '?Insulation', 'FontSize', 12, 'fontname', 'Arial', 'Location', 'southeast');
    ylabel('Insulation Score', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold');
    xlabel('Region', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold');%xL = xlim;
    colorbar
    set(gca, 'colormap', StrengthColor)
    %colormap(co)
    caxis([0 1.0])
    %colorbar('Ticks',[-5,-2,1,4,7],...
    %'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})
    %line(xL,[4 4],'Color','r');
    %save data/fig
    %exportgraphics(gcf, ['res' num2str(Resolution) '/normmean_cluster_insulation' num2str(p) '.pdf'],'ContentType','vector','BackgroundColor','none');
    %saveas(gcf, ['res' num2str(Resolution)  '/normmean_cluster_insulation' num2str(p)]);
    %get(asxes,'Position')
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'FontName', 'Arial');
    grid on
    hold off

    set(gcf, 'Position', [161 53 554 891])

    filename = [dirname '/Cluster_' num2str(k - 1) '_Insulation_' num2str(InsBoxSize) '_' ...
        'Delta_' num2str(DeltaBoxSize) '_Profiles_single_TH_' num2str(Boundary_TH) '_Reg_' ...
        num2str(RegionStart) 'to' num2str(TotalNumTADs) '.fig'];
    if s == 1
        saveas(gcf, filename)
    end
end


function u = draw_figsNorm(Mean, axValues, k, Mean_adjust, RegionStart, ...
                                TotalNumTADs,numRegions, Boundary_TH, ...
                                InsProfRange, InsulationProfNorm, ...
                                DeltaRange, DeltaInsulation, ... % DomainCenter, ...
                                Boundaries, StrengthColor, dirname, ...
                                InsBoxSize, DeltaBoxSize, s, numskip, numskipped,ChrClus)
    
    figure
    u = tiledlayout(5, 1, 'TileSpacing', 'compact');
    
    nexttile(1, [3 1])

    imagesc(Mean_adjust);
    set(gca, 'colormap', hot)
    caxis([0.6, 1.4]);
    % Uncomment the following line to preserve the X-limits of the axes
    xlim([0.5 numRegions + 0.5]); %28.5
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim([0.5 numRegions + 0.5]);
    yticks(1:5:numRegions + 1);
    yticklabels({axValues(RegionStart - numskipped:5:TotalNumTADs - numskip - numskipped)})
    xticks(1:5:numRegions + 1);
    xticklabels({axValues(RegionStart - numskipped:5:TotalNumTADs - numskip - numskipped)})
    xtickangle(90);
    colorbar;
    axis square
    RedBlue
    title(['Normalized distance Cluster ' num2str(k - 1) ' n=' num2str(length(ChrClus))]);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'FontName', 'Arial');

    nexttile(4, [2 1])
    %nexttile(10,[2 1])
    hold on
    %     ydata = RegionStart:TotalNumTADs;
    %h = yline(0, 'k-', 'LineWidth', 1);
    %         area(InsulationProfilesNorm(p+1, :), 'LineWidth', 1, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7]);
    %         r = plot(InsulationProfilesNorm(p+1, :), 'LineWidth', 1, 'Color', [0 0 0]);
    %         rr = plot(DeltaInsulation(p+1, :), 'LineWidth', 1, 'Color', [1 0 0]);
    area(InsProfRange, InsulationProfNorm(:), 'LineWidth', 1, ...
        'FaceColor', [0.7 0.7 0.7], 'EdgeColor', [0.7 0.7 0.7]);
    Ins = plot(InsProfRange, InsulationProfNorm(:), 'LineWidth', ...
        1, 'Color', [0 0 0]);

    rDelta = plot(DeltaRange, DeltaInsulation(:), 'LineWidth', 1, ...
            'Color', [1 0 0]);r = plot(DeltaRange, DeltaInsulation(:), 'LineWidth', 1, ...
         'Color', [1 0 0]);
%    plot(diff(DeltaIsulation));
    xlim([0.5 numRegions + 0.5]);
    xticks(1:5:numRegions);
    xticklabels({axValues(RegionStart - numskipped:5:TotalNumTADs - numskip - numskipped)});
    % xticklabels({round(DomainCenter(1,RegionStart:TotalNumTADs),1)});
    xtickangle(90)
    title('Insulation(Black) and delta(Red) profile');
    %grid on
    ylabel('score');
    xlabel('region');
    ylim([-1.1 1.1]);
    % yL = ylim;
    if ~isempty(Boundaries) %if there are not boundaries this won't run
        for i = 1:length(Boundaries(:, 1))
            if Boundaries(i, 2) > 1.0
                %line([Boundaries(i,1) Boundaries(i,1)],[0.52 0.55], yL, 'color',[0.6 0.6 0.6], 'linewidth', 3 );
%                 rrr = line([Boundaries(i, 1) Boundaries(i, 1)], [1.00 1.03], yL, 'color', StrengthColor(round(Boundaries(i, 2) * 100), :, :), 'linewidth', 3);
                bound = scatter(Boundaries(i, 1), 1.0, 50, StrengthColor(round(Boundaries(i, 2) * 100), :), 'v', 'filled');
            elseif Boundaries(i, 2) > Boundary_TH % this is 0.1 th of crane.
%                 rrr = line([Boundaries(i, 1) Boundaries(i, 1)], [Boundaries(i, 2) (Boundaries(i, 2) + 0.03)], yL, 'color', StrengthColor(round(Boundaries(i, 2) * 100), :), 'linewidth', 6);
                bound = scatter(Boundaries(i, 1), Boundaries(i, 2), [], StrengthColor(round(Boundaries(i, 2) * 100), :), 'v', 'filled');
            end%x-axis
        end%x-axis
        legend([Ins, bound], 'Insulation', 'Boundary', 'FontSize', 12, 'fontname', 'Arial', 'Location', 'southeast');
    else
        legend(Ins, 'Insulation', 'FontSize', 12, 'fontname', 'Arial', 'Location', 'southeast');
    end
    
    title('Insulation and delta profile', 'FontSize', 12, 'fontname', 'Arial', 'VerticalAlignment', 'bottom');
%     legend([r, rr], 'Insulation', '?Insulation', 'FontSize', 12, 'fontname', 'Arial', 'Location', 'southeast');
    ylabel('Insulation Score', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold');
    xlabel('Region', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold');%xL = xlim;
    colorbar
    set(gca, 'colormap', StrengthColor)
    %colormap(co)
    caxis([0.1 1.0])
    %colorbar('Ticks',[-5,-2,1,4,7],...
    %'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})
    %line(xL,[4 4],'Color','r');
    %save data/fig
    %exportgraphics(gcf, ['res' num2str(Resolution) '/normmean_cluster_insulation' num2str(p) '.pdf'],'ContentType','vector','BackgroundColor','none');
    %saveas(gcf, ['res' num2str(Resolution)  '/normmean_cluster_insulation' num2str(p)]);
    %get(asxes,'Position')
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'FontName', 'Arial');
    grid on
    hold off

    set(gcf, 'Position', [161 53 554 891])

    filename = [dirname '/Norm_Cluster_' num2str(k - 1) '_Insulation_' num2str(InsBoxSize) '_' ...
        'Delta_' num2str(DeltaBoxSize) '_Profiles_single_TH_' num2str(Boundary_TH) '_Reg_' ...
        num2str(RegionStart) 'to' num2str(TotalNumTADs) '.fig'];
    if s == 1
        saveas(gcf, filename)
    end
end

