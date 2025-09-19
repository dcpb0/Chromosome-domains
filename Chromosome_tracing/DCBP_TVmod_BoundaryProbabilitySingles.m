%% Single Trace Boundary Probability Calculation – Nov 2024
% Author(s): DCPB (original), improved by T. Verheijen.
% Description:
% This script calculates the proportions of boundary locations in single-cell 
% traces. Insulation is computed following Crane et al. (2015) using square 
% boxes along the diagonal. Boundary prevalence is estimated similarly to 
% Chen et al. (2021): data are first linearly filled, then insulation is 
% computed (as in Crane), and finally the number of detected boundaries is 
% normalized by the total number of single traces.
%
% Notes:
% The insulation and delta box size can be adjusted to detect more or 
% fewer changes.
% The threshold for boundary detection should be set carefully:
%   • Crane et al. (2015) used 0.1, which works well for averaged data.  
%   • For single traces, this threshold is too low; a minimum of 0.5 is 
%     recommended (tested on chrV).

close all
clearvars -except Chr y_fit

addpath('/scicore/home/mangos/pulido0000/Tracing/BasicFunctions')

% DECLARE if you would like to save the figures by setting savefigs = 1,
% else 0
% Note that this will also save all figures of the singles (many figures)
savefigs = 1;

% load ax values for axis labels
load('1to40_ce10_axValues.mat')
axValues  = axValues(14:40);

%Input this if you want to normalize the insulation to the average
% MeanInsulationAverage=  0.9515; %WT
%MeanInsulationAverage=  0.9486; %MT

% load the Y_Fit of the average
%load('R15_N2_2to40_YFit_min14.mat');

% load Chr data
%load('old_data/AllAges_AllChromosomes_merge_N2vsMET2.mat')

% Create colormap for Boundary Strength
StrengthColor = flipud(parula(450));
StrengthColor = StrengthColor(80:430, :);

% Define constants:
RegionStart = 1;
minRegions = 25;
Boundary_TH = 0.5;
TotalNumTADs = 43-13;
SkipTADs = sort([28, 29, 40], 'descend') -13;
SkipTADs = SkipTADs(SkipTADs < TotalNumTADs);
numRegions = TotalNumTADs - RegionStart - length(SkipTADs) + 1;

% Don't forget to specify insulation and delta sizes
InsBoxSize = 3;
DeltaBoxSize = 2;

% silence if you don't need to select an age range
AgeStart = 2; %change these to specify embryo age
AgeStop = 140;
AgeRange = ([ num2str(AgeStart) 'to' num2str(AgeStop) 'cell']); 

n=0;
for k = 1:length(Chr)
    if  Chr(k).Age >= AgeStart && Chr(k).Age <= AgeStop
        n=n+1;
        Chr2(n)= Chr(k);
    end
end
Chr = Chr2;

%
name = 'Insulation_singles_2to140';
dirname = ([num2str(name) '_min_' num2str(minRegions) '_TH_' ...
    num2str(Boundary_TH) '_Reg_' num2str(RegionStart) 'to' ...
    num2str(TotalNumTADs)]);
if ~isfolder(dirname) && savefigs == 1
    mkdir(dirname);
    mkdir([dirname '/singles'])
end

% Select quality traces with high region counts
Chr_chosen = Chr(arrayfun(@(x) sum(x.r), Chr) >= minRegions);
Chr = Chr_chosen;

% saveas(gcf, [dirname '/2to40_Num_Regions_Distribution_singles.fig'])

% Used to normalize matrices with y_fit of the average
% pre compute the genomic distances for normalization

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

% These two ranges are the x-indices which correspond to the correct
% regions. If you don't use these ranges, all of you information will
% appear shifted and not make sense.
InsProfRange = InsBoxSize + 1:DomainSize - InsBoxSize;
DeltaRange = InsBoxSize + DeltaBoxSize + 1:DomainSize - ...
                                                InsBoxSize - DeltaBoxSize;
% Preallocate space
AllStrength = cell(length(Chr_chosen), 1);
AllBoundaries = cell(length(Chr_chosen), 1);
tracesmissing = zeros(length(Chr_chosen), 1) - 1;

numskipped = 0;
while sum(RegionStart == SkipTADs) ~= 0
    SkipTADs(end) = [];
    RegionStart = RegionStart + 1;
    numskipped = numskipped + 1;
end

numRegions = TotalNumTADs - RegionStart + 1;

% Find the mean of all chosen traces to fill missing values

MeanofMeans = zeros(TotalNumTADs);

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
            %                 if dist <= 6 && ~isnan(dist)
            %                     DisList(k) = dist;
            %                 end
        end
        DisList(DisList == -1) = [];

        MeanofMeans(i, j) = mean(DisList);
    end
end

for i = 1:length(SkipTADs)
    MeanofMeans(SkipTADs(i), :) = [];
    MeanofMeans(:, SkipTADs(i)) = [];
end

%% Main Loop: Compute Boundary positions and strength
for k = 1:length(Chr_chosen)
    numRegions = TotalNumTADs - RegionStart + 1;

    Mean = zeros(numRegions);

    for i = 1:numRegions
        for j = 1:numRegions
            if Chr_chosen(k).r(i) == 1 && Chr_chosen(k).r(j) == 1
                %Dist(i,j) = distance(C, i, j);
                xsq = (Chr_chosen(k).x(i) - Chr_chosen(k).x(j))^2;
                ysq = (Chr_chosen(k).y(i) - Chr_chosen(k).y(j))^2;
                zsq = (Chr_chosen(k).z(i) - Chr_chosen(k).z(j))^2;
                Mean(i,j) = sqrt(xsq + ysq + zsq);
            else
                Mean(i,j) = NaN;
            end
        end
    end
    for i = 1:length(SkipTADs)
        Mean(SkipTADs(i), :) = [];
        Mean(:, SkipTADs(i)) = [];
    end

    numRegions = length(Mean);

    nondetected = zeros(numRegions - 2, 1);
    for i = 1:numRegions - 2
        % Check for three or more consecutive missing columns
        nondetected(i) = isnan(mean(Mean(i:i + 2, :), 'all', 'omitnan'));
    end

    % If there are any sections of three or more missing columns, we move
    % on to the next trace
    if sum(nondetected) ~= 0
        continue;
    end

    % Check for first and last rows containing missing values (can lead to
    % strange linear filling if left untouched
    if isnan(mean(Mean(1, :), 'all', 'omitnan'))
        Mean(1, :) = MeanofMeans(1, :);
        Mean(:, 1) = MeanofMeans(:, 1);
    end

    if isnan(mean(Mean(end, :), 'all', 'omitnan'))
        Mean(end, :) = MeanofMeans(end, :);
        Mean(:, end) = MeanofMeans(:, end);
    end

    % Fill the Mean matrix with Linear Interpolation
    Mean_filtered = fillmissing(Mean, 'linear'); % Fill Rows
    Mean_filtered = fillmissing(Mean_filtered, 'linear', 2); % Fill Columns

    negatives = find(Mean_filtered < 0);
    if ~isempty(negatives)
        Mean_filtered(negatives) = MeanofMeans(negatives);
    end

    % If the Genomic Normalization is being calculated during each loop, 
    % uncomment this section and comment out the Mean_adjust loop below
    %
    %     Dis_pxl = Mean_filtered(tril(true(size(Mean_filtered)), -1));
    %     R = corrcoef(GenomicDis_pxl, Dis_pxl);
    %     x = log(GenomicDis_pxl);
    %     y = log(Dis_pxl);
    %     X = [ones(size(x)), x];
    %     [b, bint] = regress(y, X);
    %     x = GenomicDis_pxl';
    %     y_fit = exp(b(1)) * x.^b(2);
    %     norms = Dis_pxl./y_fit';
    %     Mean_adjust = squareform(norms) + eye(numRegions);


    % Perform Genomic Normalization based with read in y_fit. Comment this
    % loop out if doing individual Genomic Normalization

    Mean_adjust = Mean_filtered;
    Dis_pxl = Mean_filtered(tril(true(size(Mean_filtered)), -1));
    norms = Dis_pxl./y_fit';
    Mean_adjust = squareform(norms) + eye(numRegions);

    % Use this if the y_fit is done in the whole matrix. No tril!
    
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
    InsProfNorm = log2(InsProfile/mean(InsProfile));

    %Normalize the insulation to the average (Silence the one before)
   %InsProfNorm = log2(InsProfile/MeanInsulationAverage);

    % Find the change in insulation values
    DeltaIns = zeros(length(InsProfNorm) - 2 * DeltaBoxSize, 1);
    % We include -2 * DeltaBoxSize, because we cannot use the two
    % DeltaBoxSize means on either end of the InsProfNorm array
    for i = 1:length(DeltaIns)
        DeltaIns(i) = mean(InsProfNorm(i:i + DeltaBoxSize - 1)) - ...
            mean(InsProfNorm(i + DeltaBoxSize + 1:i + 2 * DeltaBoxSize));
    end

    % Find peaks and valleys of the insulation
    [InsMax, InsMaxIdx] = findpeaks(InsProfNorm);
    [InsMin, InsMinIdx] = findpeaks(-InsProfNorm);
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
    if isempty(Boundaries)
        continue;
    end

    % Filter for boundaries above the boundary threshold
    Boundaries = Boundaries(Boundaries(:, 2) > Boundary_TH, :);

    [~, ind, ~] = unique(Boundaries(:, 2));
    Boundaries = Boundaries(sort(ind, 'ascend'), :); % Ensure lo-high order

    AllBoundaries{k} =  Boundaries(:, 1);
    AllStrength{k} =  Boundaries(:, 2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Draw figures for this trace
%     t = draw_figs(Mean, axValues, k, Mean_filtered, ...
%                                 RegionStart, numRegions, Boundary_TH, ...
%                                 InsProfRange, InsProfNorm, ...
%                                 DeltaRange, DeltaIns, ...% DomainCenter, 
%                                 Boundaries, StrengthColor, dirname, ...
%                                 InsBoxSize, DeltaBoxSize);

    t = draw_figs(Mean, axValues, k, Mean_filtered, ...
        RegionStart, TotalNumTADs, numRegions, Boundary_TH, ...
        InsProfRange, InsProfNorm, ...
        DeltaRange, DeltaIns, ...% DomainCenter,
        Boundaries, StrengthColor, dirname, InsBoxSize, DeltaBoxSize, ...
        savefigs, length(SkipTADs), numskipped);
    
    close all
end % End for 1:length(Chr_chosen)

AllBoundaries = AllBoundaries(~cellfun('isempty', AllBoundaries));
AllStrength = AllStrength(~cellfun('isempty', AllStrength));

total_occurrences = zeros(numRegions, 1);
for i = 1:numRegions
    Counts = cellfun(@(x) find(x == i), AllBoundaries, 'UniformOutput', 0);
    total_occurrences(i) = numel([Counts{:}]);
end

% The new line is uncommented above the old line
total_occurrences = total_occurrences/length(AllBoundaries);

if savefigs == 1
    save([dirname '/AllBoundaries_TH_' num2str(Boundary_TH) '_min_' num2str(minRegions) '.mat'], 'AllBoundaries');
    save([dirname '/AllStrength_TH_' num2str(Boundary_TH) '_min_' num2str(minRegions) '.mat'], 'AllStrength');
    save([dirname '/Probability_TH_' num2str(Boundary_TH) '_min_' num2str(minRegions) '.mat'], 'total_occurrences');
end


figure
hold on
ydata = 1:numRegions; %this is the x data
sz = numRegions;
r = scatter(ydata, total_occurrences, sz, 'filled', 'MarkerFaceColor', [0.6350 0.0780 0.1840], 'MarkerFaceAlpha', 0.8);
rr= plot(total_occurrences(:), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
%rr= plot(total_occurrences(:), 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
yline(mean(nonzeros(total_occurrences)), '--', 'Mean', 'LineWidth', 3);
grid on
yline(mean(nonzeros(total_occurrences)) + std(total_occurrences), ':', 'LineWidth', 2);
yline(mean(nonzeros(total_occurrences)) - std(total_occurrences), ':', 'LineWidth', 2);
xlim([1 numRegions]);
xticks(1:5:numRegions);
xticklabels({axValues(RegionStart - numskipped:5:TotalNumTADs - length(SkipTADs) - numskipped)});
%xticklabels({round(DomainCenter(1,RegionStart:TotalNumTADs),1)})
xtickangle(90);
ylim([0 0.25])
title('Boundary probability from single cells ', 'FontSize', 12, 'fontname', 'Arial', 'VerticalAlignment', 'bottom');
ylabel('Boundary probability', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold');
xlabel('Region', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold')%xL = xlim;
set(gca, 'LineWidth', 1.5, 'FontSize', 12, 'FontName', 'Arial');
hold off

% Remove zeros and NaN values
filtered_data = total_occurrences(total_occurrences ~= 0 & ~isnan(total_occurrences));
% Calculate the mean
mean_value = mean(filtered_data);
% Plot the mean as a horizontal line
yline(mean_value, '--', 'Mean', 'LineWidth', 3);

if savefigs == 1
    saveas(gcf, [dirname '/Boundary_Probability_Insulation_' ...
        num2str(InsBoxSize) '_Delta_' num2str(DeltaBoxSize) ...
        '_Reg_' num2str(RegionStart) 'to' num2str(TotalNumTADs) '.fig'])
end

%
numericArray = cell2mat(AllStrength);

figure
numBins = 10;
histogram(numericArray, 'Normalization', 'probability','NumBins', numBins, 'FaceColor', [0.1660 0.5740 0.0880]);
title('Boundary Probability vs Strength');
xlabel('Boundary Strength');
ylabel('Probability');
saveas(gcf, [dirname '/ProbabilyVSstrength_' ...
        num2str(InsBoxSize) '_Delta_' num2str(DeltaBoxSize) ...
        '_Reg_' num2str(RegionStart) 'to' num2str(TotalNumTADs) '.fig'])


%% Declare Functions

% Mean, axValues, k, Mean_filtered, RegionStart, ...
%                                 numRegions, Boundary_TH, ...
%                                 InsProfRange, InsulationProfNorm, ...
%                                 DeltaRange, DeltaInsulation, ... % DomainCenter, ...
%                                 Boundaries, StrengthColor, dirname, ...
%                                 InsBoxSize, DeltaBoxSize ...

function t = draw_figs(Mean, axValues, k, Mean_filtered, RegionStart, ...
                                TotalNumTADs,numRegions, Boundary_TH, ...
                                InsProfRange, InsulationProfNorm, ...
                                DeltaRange, DeltaInsulation, ... % DomainCenter, ...
                                Boundaries, StrengthColor, dirname, ...
                                InsBoxSize, DeltaBoxSize, s, numskip, numskipped)
    figure
    t = tiledlayout(5, 2, 'TileSpacing', 'compact');

    nexttile(1, [3 1])
    imagesc(Mean);
    set(gca, 'colormap', hot)
    caxis([0, 1.3])
    % Uncomment the following line to preserve the X-limits of the axes
    xlim([0.5  numRegions + 0.5]); %28.5
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim([0.5  numRegions + 0.5]);
    yticks(1:5:numRegions + 1)
    yticklabels({axValues(RegionStart - numskipped:5:TotalNumTADs - numskip - numskipped)})
    xticks(1:5:numRegions + 1)
    xticklabels({axValues(RegionStart - numskipped:5:TotalNumTADs - numskip - numskipped)})
    xtickangle(90)
    colorbar;
    RedBlue
    title(['Distance matrix trace: ' num2str(k)])
    axis square
    %set(gcf,'units','centimeters','position',[10,10,30,30])
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'FontName', 'Arial');

    nexttile(2, [3 1])
    imagesc(Mean_filtered);
    set(gca, 'colormap', hot)
    caxis([0, 1.3]);
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
    title(['Filled distance trace: ' num2str(k)]);
    set(gca, 'LineWidth', 1.5, 'FontSize', 10, 'FontName', 'Arial');

    nexttile(8, [2 1])
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
    Delta = plot(DeltaRange, DeltaInsulation(:), 'LineWidth', 1, ...
        'Color', [1 0 0]);
    %plot(diff(DeltaIsulation));
    xlim([0.5 numRegions + 0.5]);
    xticks(1:2:numRegions + 1);
    xticklabels({axValues(RegionStart - numskipped:2:TotalNumTADs - numskip - numskipped)});
    % xticklabels({round(DomainCenter(1,RegionStart:TotalNumTADs),1)});
    xtickangle(90)
    title('Insulation(Black) and delta(Red) profile');
    %grid on
    ylabel('score');
    xlabel('region');
    ylim([-3.5 3.5]);
    % yL = ylim;
    if ~isempty(Boundaries) %if there are not boundaries this won't run
        for i = 1:length(Boundaries(:, 1))
            if Boundaries(i, 2) > 3.5
                %line([Boundaries(i,1) Boundaries(i,1)],[0.52 0.55], yL, 'color',[0.6 0.6 0.6], 'linewidth', 3 );
                bound = scatter(Boundaries(i, 1), 3.47, 50, StrengthColor(350, :), 'v', 'filled');
            elseif Boundaries(i, 2) > Boundary_TH % this is 0.1 th of crane.
                bound = scatter(Boundaries(i, 1), Boundaries(i, 2), [], StrengthColor(round(Boundaries(i, 2) * 100), :), 'v', 'filled');
            end%x-axis
        end%x-axis
        legend([Ins, Delta, bound], 'Insulation', 'Delta', 'Boundary', 'FontSize', 12, 'fontname', 'Arial', 'Location', 'southeast');
    else
        legend([Ins, Delta], 'Insulation', 'Delta', 'FontSize', 12, 'fontname', 'Arial', 'Location', 'southeast');
    end
    title('Insulation and delta profile', 'FontSize', 12, 'fontname', 'Arial', 'VerticalAlignment', 'bottom');
    ylabel('Score', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold');
    xlabel('Region', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold');%xL = xlim;
    colorbar
    set(gca, 'colormap', StrengthColor)
    %colormap(co)
    caxis([0.1 3.5])
    %colorbar('Ticks',[-5,-2,1,4,7],...
    %'TickLabels',{'Cold','Cool','Neutral','Warm','Hot'})
    %line(xL,[4 4],'Color','r');
    %save data/fig
    %exportgraphics(gcf, ['res' num2str(Resolution) '/normmean_cluster_insulation' num2str(p) '.pdf'],'ContentType','vector','BackgroundColor','none');
    %get(asxes,'Position')
    set(gca, 'FontName', 'Arial', 'FontSize', 10);
    grid on
    hold off

    filename = [dirname '/singles/Singles_' num2str(k) '_Insulation_' num2str(InsBoxSize) '_' ...
        'Delta_' num2str(DeltaBoxSize) '_Profiles_single_TH_' num2str(Boundary_TH) '_Reg_' ...
        num2str(RegionStart) 'to' num2str(TotalNumTADs) '.fig'];
    if s == 1
        saveas(gcf, filename)
    end
end



