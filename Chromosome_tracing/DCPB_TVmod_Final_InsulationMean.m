%% Single Trace Boundary Probability Calculation 
% Author(s): DCPB (original), improved by T. Verheijen
%
% Insulation is calculated following Crane et al. (2015) using square boxes
% along the diagonal 


close all
clearvars -except Chr

addpath('/scicore/home/mangos/pulido0000/Tracing/BasicFunctions')

% load ax values for axis labels
load('1to40_ce10_axValues.mat')

% load the Y_Fit
% load('old_data/Min14_2to40cells_allchr_Y_Fit.mat');

% load the Chr data
% load('old_data/2_40_AllChromosomesNoNucL_merge4_15_17.mat')

ColorMap = load('RedBlue.txt');
redblue = ColorMap/255;

% Create colormap for Boundary Strength
StrengthColor = flipud(parula(450));
StrengthColor = StrengthColor(80:430, :);

% Define constants:
RegionStart = 1;
minRegions = 1;
Boundary_TH = 0.1;
TotalNumTADs = 13;
SkipTADs = sort([], 'descend'); %28, 29, 40
numRegions = TotalNumTADs - length(SkipTADs);
InsBoxSize = 2;
DeltaBoxSize = 1;

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


%%

name = 'MT_insulation_median_R15_2to140_';
dirname = ([num2str(name) '_min_' num2str(minRegions) '_TH_' ...
    num2str(Boundary_TH) '_TadNum_' num2str(numRegions)]);
if ~isfolder(dirname)
    mkdir(dirname);
end

% Select quality Cells with high region counts
Chr_chosen = Chr(arrayfun(@(x) sum(x.r), Chr) >= minRegions);
Chr = Chr_chosen;

figure
histogram(arrayfun(@(x) sum(x.r), Chr_chosen))
title('Number of Regions Distribution')
xlabel('Number of regions detected')
ylabel('Number of traces')

filename = [dirname '/' num2str(RegionStart) 'to' num2str(TotalNumTADs)...
    '_Num_Regions_Distribution_Mean_min_' num2str(minRegions) '.fig'];
if ~isfile(filename)
    saveas(gcf, filename);
end

% Used to normalize matrices with y_fit of the average

fid = fopen('Probes_startends2829_40_5s_rDNA_Ce10.txt', 'r');
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

InsProfRange = InsBoxSize + 1:DomainSize - InsBoxSize;
DeltaRange = InsBoxSize + DeltaBoxSize + 1:DomainSize - ...
InsBoxSize - DeltaBoxSize;

%%
% Find the mean of all chosen traces
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
%                 if dist <= 6 && ~isnan(dist)
%                     DisList(k) = dist;
%                 end
            end
            DisList(DisList == -1) = [];

            Mean(i, j) = median(DisList); %ATTENTION:i changed for the median 
        end
    end

    for i = 1:length(SkipTADs)
        Mean(SkipTADs(i), :) = [];
        Mean(:, SkipTADs(i)) = [];
    end

%% Main Section: Compute Boundary positions and strength

% Fill the Mean matrix with Linear Interpolation
Mean_filtered = fillmissing(Mean, 'linear'); % Fill Rows
Mean_filtered = fillmissing(Mean_filtered, 'linear', 2); % Fill Columns

negatives = find(Mean_filtered < 0);
if ~isempty(negatives)
    Mean_filtered(negatives) = mean(Mean_filtered, 'all', 'omitnan');
end

% If the Genomic Normalization is being calculated during each loop,
% uncomment this section and comment out the Mean_adjust loop below
%
Dis_pxl = Mean_filtered(tril(true(size(Mean_filtered)), -1));
R = corrcoef(GenomicDis_pxl, Dis_pxl);
x = log(GenomicDis_pxl);
y = log(Dis_pxl);
X = [ones(size(x)), x];
[b, bint] = regress(y, X);
x = GenomicDis_pxl';
y_fit = exp(b(1)) * x.^b(2);
norms = Dis_pxl./y_fit';
Mean_adjust = squareform(norms) + eye(numRegions);

% Perform Genomic Normalization based with read in y_fit. Comment this
% loop out if doing individual Genomic Normalization
%     Mean_adjust = Mean_filtered;
%     n = 1;
%     for i = 1:TotalNumTADs
%         for j = 1:TotalNumTADs
%             if i~=j
%                 Mean_adjust(i,j) = Mean_filtered(i,j)/y_fit(n);
%                 n = n + 1;
%             else
%                 Mean_adjust(i,j) = 1;
%             end
%         end
%     end

% figure
% imagesc(Mean_adjust);
% set(gca, 'colormap', redblue)
% caxis([0.5 1.5])

%%

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

% If no valid boundaries, we simply skip to making the figures with no
% boundaries.
if ~isempty(Boundaries)
    % Filter for boundaries above the boundary threshold
    Boundaries = Boundaries(Boundaries(:, 2) > Boundary_TH, :);

    [~, ind, ~] = unique(Boundaries(:, 2));
    Boundaries = Boundaries(sort(ind, 'ascend'), :); % Ensure lo-high order
end
%%
save([dirname '/InsulationProfileNorm_' ...
    num2str(InsBoxSize) '_Delta_' num2str(DeltaBoxSize) '_'...
    'Profile_TH_' num2str(Boundary_TH) '_Reg_' ...
    num2str(RegionStart) 'to' num2str(TotalNumTADs) '.mat'], 'InsProfNorm');

save([dirname '/Boundaries_Ins_' ...
    num2str(InsBoxSize) '_Delta_' num2str(DeltaBoxSize) '_'...
    'Profile_TH_' num2str(Boundary_TH) '_Reg_' ...
    num2str(RegionStart) 'to' num2str(TotalNumTADs) '.mat'], 'Boundaries');

save([dirname '/AllChr_' num2str(AgeStart) 'to' num2str(AgeStop) 'cells.mat'], 'Chr');
save([dirname '/yFit_' num2str(AgeStart) 'to' num2str(AgeStop) 'cells.mat'], 'y_fit');

%% Draw figures for the mean matrix and boundaries
t = draw_figs(Mean, axValues, Mean_filtered, RegionStart, numRegions,...
        Boundary_TH, InsProfRange, InsProfNorm, DeltaRange, DeltaIns, ...% DomainCenter, ...
        Boundaries, StrengthColor, dirname, InsBoxSize, DeltaBoxSize);


%% Declare Functions

% function [DS, DE, DC] = read_domains(TotalNumTADs)
%     fid = fopen('~/BasicFunctions/Probes_startends2829_40_5s_rDNA_Ce10.txt', 'r');
%     tline = fgetl(fid);
%     linenum = 1;
%     DS = zeros(100, 1) - 1;
%     DE = zeros(100, 1) - 1;
%     DC = zeros(100, 1) - 1;
%     while ischar(tline)
%         C = strsplit(tline);
%         DS(linenum) = str2double(C{2});
%         DE(linenum) = str2double(C{3});
%         DC(linenum) = (DS(linenum) + DE(linenum))/2.0;
%         tline = fgetl(fid);
%         linenum = linenum + 1;
%     end
%     fclose(fid);
%     
%     DS = DS(1:TotalNumTADs);
%     DE = DE(1:TotalNumTADs);
%     DC = DC(1:TotalNumTADs);
% 
%     DC = DC./1000000; %convert from bp to Mb
% end
% 
% function generate_region_histogram(a)
%     figure
%     histogram(arrayfun(@(x) sum(x.r), a))
%     title('Number of Regions Distribution')
%     xlabel('Number of regions detected')
%     ylabel('Number of traces')
%     ColorMap = load('~/BasicFunctions/RedBlue.txt');
%     colormap(ColorMap/255);
% end
% 
% function Mean = computeOverallMean(C, SkipTADs, NumRegions)
%     Mean = zeros(NumRegions);
%     
%     for i = 1:NumRegions
%         for j = 1:NumRegions
%             DisList = zeros(1, length(C)) - 1;
%             for k = 1:length(C)
%                 if ~C(k).r(i) || ~C(k).r(j)
%                     continue;
%                 end
%     
%                 dist = distance(C(k), i, j);
%                 % Confirm distance ≤ 6 and ensures we don't have NaN's
%                 if dist <= 6 && ~isnan(dist)
%                     DisList(k) = dist;
%                 end
%             end
%             DisList(DisList == -1) = [];
%     
%             Mean(i, j) = mean(DisList);
%         end
%     end
%     
%     Mean = remove_TADs(Mean, SkipTADs);
% end
% 
% function Dist = computeDistMat(C, SkipTADs, NumRegions)
%     Dist = zeros(NumRegions);
% 
%     for i = 1:NumRegions
%         for j = 1:NumRegions
%             if C.r(i) == 1 && C.r(j) == 1
%                 Dist(i,j) = distance(C, i, j);
%             else
%                 Dist(i,j) = NaN;
%             end
%         end
%     end
% 
%     Dist = remove_TADs(Dist, SkipTADs);
% end
% 
% function mat = remove_TADs(mat, skips)
%     for i = 1:length(skips)
%         mat(skips(i), :) = [];
%         mat(:, skips(i)) = [];
%     end
% end
% 
% function d = distance(vect, i, j)
% % Compute distance given a 3-dimensional vector as input
% xsq = (vect.x(i) - vect.x(j))^2;
% ysq = (vect.y(i) - vect.y(j))^2;
% zsq = (vect.z(i) - vect.z(j))^2;
% d = sqrt(xsq + ysq + zsq);
% end

function t = draw_figs(Mean, axValues, Mean_filtered, RegionStart, ...
                                TotalNumTADs, Boundary_TH, ...
                                InsProfRange, InsulationProfNorm, ...
                                DeltaRange, DeltaInsulation, ... % DomainCenter, ...
                                Boundaries, StrengthColor, dirname, ...
                                InsBoxSize, DeltaBoxSize)
    
        figure
        t = tiledlayout(5, 1, 'TileSpacing', 'compact');
    
        nexttile(1, [3 1])
        imagesc(Mean_filtered);
        set(gca, 'colormap', hot)
        RedBlue
        caxis([0, 1.5])
        % Uncomment the following line to preserve the X-limits of the axes
        xlim([RegionStart - 0.5 TotalNumTADs + 0.5]); %28.5
        % Uncomment the following line to preserve the Y-limits of the axes
        ylim([RegionStart - 0.5 TotalNumTADs + 0.5]);
        yticks(RegionStart:2:TotalNumTADs)
        yticklabels({axValues(RegionStart:2:TotalNumTADs)})
        xticks(RegionStart:2:TotalNumTADs)
        xticklabels({axValues(RegionStart:2:TotalNumTADs)})
        xtickangle(90)
        colorbar;
        title(['Mean of Distance Matrices - Regions ' ...
                    num2str(RegionStart) 'to' num2str(TotalNumTADs)])
        axis square
        %set(gcf,'units','centimeters','position',[10,10,30,30])
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
        Delta = plot(DeltaRange, DeltaInsulation(:), 'LineWidth', 1, ...
            'Color', [1 0 0]);
        %plot(diff(DeltaIsulation));
        xlim([RegionStart-0.5 TotalNumTADs+0.5]);
        xticks(1:2:TotalNumTADs);
        xticklabels({axValues(RegionStart:2:TotalNumTADs)});
        % xticklabels({round(DomainCenter(1,RegionStart:TotalNumTADs),1)});
        xtickangle(90)
        title('Insulation and delta profile');
        %grid on
        ylabel('score');
        xlabel('region');
        ylim([-0.55 0.55]);
        % yL = ylim;
        if ~isempty(Boundaries) %if there are not boundaries this won't run
            for i = 1:length(Boundaries(:, 1))
                if Boundaries(i, 2) > 1.0
                    %line([Boundaries(i,1) Boundaries(i,1)],[0.52 0.55], yL, 'color',[0.6 0.6 0.6], 'linewidth', 3 );
                    bound = scatter(Boundaries(i, 1), 1.0, 50, StrengthColor(round(Boundaries(i, 2) * 100), :), 'v', 'filled');
                elseif Boundaries(i, 2) > Boundary_TH % this is 0.1 th of crane.
                    bound = scatter(Boundaries(i, 1), Boundaries(i, 2), 50, StrengthColor(round(Boundaries(i, 2) * 100), :), 'v', 'filled');
                end%x-axis
            end%x-axis
            legend([Ins, Delta, bound], 'Insulation', 'Delta', 'Boundary', 'FontSize', 12, 'fontname', 'Arial', 'Location', 'southeast');
        else
            legend([Ins, Delta], 'Insulation', 'Delta', 'FontSize', 12, 'fontname', 'Arial', 'Location', 'southeast');
        end

        title('Insulation and delta profile', 'FontSize', 12, 'fontname', 'Arial', 'VerticalAlignment', 'bottom');
        ylabel('Insulation Score', 'FontSize', 12, 'fontname', 'Arial', 'fontweight', 'bold');
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
        %saveas(gcf, ['res' num2str(Resolution)  '/normmean_cluster_insulation' num2str(p)]);
        %get(asxes,'Position')
        set(gca, 'FontName', 'Arial', 'FontSize', 10);
        grid on
        hold off

        set(gcf, 'Position', [161 53 554 891])

        filename = [dirname '/mean_of_means_Insulation_' ...
            num2str(InsBoxSize) '_Delta_' num2str(DeltaBoxSize) '_'...
                'Profile_TH_' num2str(Boundary_TH) '_Reg_' ...
                num2str(RegionStart) 'to' num2str(TotalNumTADs) '.fig'];
        if ~isfile(filename)
            saveas(gcf, filename)
        end
end