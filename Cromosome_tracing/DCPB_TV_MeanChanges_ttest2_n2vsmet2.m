%% This script creates a heatmap of p-values for the Mean changes N2vsMET2
% Here we use the ttest2 function to find the statistic for the difference
% in means at each cell of the DisList variables for the N2 and MET2
% datasets. We are able to use ttest2 as it does not require the data
% inputs to be of the same dimensions and if we specify 'unequal' it does
% not assume equal variances and simply computes the Welch's t-test stat.

% This script was created by TV and DCPB 

clear
% close all

bfmatlabpath ='/scicore/home/mangos/pulido0000/Tracing/BasicFunctions';
pathString = genpath(bfmatlabpath);
addpath(pathString)

% load ax values for axis labels
load('1to40_ce10_axValues.mat')
axValues  = axValues(14:40);

% clear N2DisListMerge MET2DisListMerge

TotalTADNum = 30;

AgeStart = 2; %change these to specify embryo age
AgeStop = 140;
AgeRange = ([ num2str(AgeStart) 'to' num2str(AgeStop) 'cell']);

name = ['N2_TRYmin10_' AgeRange];
minRegions = 7; %Remove traces with less than the minRegions

filename{1,1}= '/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/01_N2/N2_R15_30REG_AllAges_merge4_15_17.mat';
filename{1,2}= '/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/02_GW_met2set25/GW_30REG_AllAges_merge10_14.mat';
filename{1,3}= '/scicore/home/mangos/pulido0000/Tracing/CondensinMutants/Allchr_Condensin/PULIDO2024/03_CondensinMT/MT_30REG_Allages_merge25_26_27.mat';

for p= 1:3;

    %clearvars -except Mean DisListAll TH AgeStart AgeStop AgeRange name minRegions N2 MT GW p filename Counts ps
    load([filename{1,p}])

    n= length(Chr);

    n=0;
    for k = 1:length(Chr)
        if  Chr(k).Age >= AgeStart && Chr(k).Age <= AgeStop
            n=n+1;
            Chr2(n)= Chr(k);
        end
    end

    Chr= Chr2;

    Chr_chosen = Chr(arrayfun(@(x) sum(x.r), Chr) >= minRegions);
    Chr = Chr_chosen;

    NumRegions = length(Chr(1).x);

    %calc and plot mean
    for i = 1:NumRegions
        for j = 1:NumRegions
            DisList = [];
            for k = 1:length(Chr)
                if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                    if ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5 <=6
                        DisList = [DisList ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5];
                    end
                end
            end
            Mean{p}(i,j) = median(DisList);
            DisListAll{p}{i,j} = DisList;
        end

    end
    SkipTADs= sort([28,29,40], 'descend')-13;

    for i =1:length(SkipTADs)
        Mean{p}([SkipTADs(i)],:)=[];
        Mean{p}(:,[SkipTADs(i)])=[];
        DisListAll{p}([SkipTADs(i)],:)=[];
        DisListAll{p}(:,[SkipTADs(i)])=[];
    end

  
end

N2DisListAll =  DisListAll{1};
MET2DisListAll= DisListAll{2};
CONDisListAll = DisListAll{3};

N2Mean =  Mean{1};
MET2Mean= Mean{2};
CONDMean = Mean{3};

Diff = MET2Mean - N2Mean;
Diff2= CONDMean - N2Mean;

% Define the RGB values for the 5 colors
cmap = [
    0, 0, 0.2;  % Blue dark
    0, 0.2, 0.2;  % Green
    0.2, 0.4, 0.3;  % Green
    1, 1, 1;  % White (Middle)
    0.2, 0.4, 0.6;  % blue
    0.2, 0.2, 0.4;  % blue
    0, 0, 0   % Black
    ];

% Create the colormap by interpolating the colors
n = 256;  % Number of colors in the final colormap
custom_cmap = interp1(linspace(0, 1, size(cmap, 1)), cmap, linspace(0, 1, n));

% Apply the colormap
colormap(custom_cmap);

figure
imagesc(Diff)
colormap(custom_cmap);
colorbar
axis square
title(['Significant mean changes N2 < MET2']);
caxis([-0.25 0.25]);
yticks(1:5:TotalTADNum)
yticklabels({axValues(1:5:TotalTADNum)})
xticks(1:5:TotalTADNum)
xticklabels({axValues(1:5:TotalTADNum)})
xtickangle(90)
%saveas(gcf,'Medians Met2 - N2','fig');
%print('Anova bin1toLate p05','-dsvg');

figure
imagesc(Diff2)
colormap(custom_cmap);
colorbar
axis square
title(['Significant mean changes N2 < Condensin I']);
caxis([-0.25 0.25]);
yticks(1:5:TotalTADNum)
yticklabels({axValues(1:5:TotalTADNum)})
xticks(1:5:TotalTADNum)
xticklabels({axValues(1:5:TotalTADNum)})
xtickangle(90)
%saveas(gcf,'Medians Met2 - N2','fig');
%print('Anova bin1toLate p05','-dsvg');

data_N2= N2Mean(N2Mean ~= 0);
data_Met = MET2Mean(MET2Mean ~= 0);
data_Con = CONDMean(CONDMean ~= 0);

%'Normalization', 'pdf': Normalizes the histogram so that it represents a
% probability density function (PDF). This makes the histogram's y-axis comparable to
% the density plot.

figure
numBins = 10;
histogram(data_N2, 'Normalization', 'pdf','FaceColor', [0 0.4470 0.7410], 'NumBins', numBins);
hold on;
histogram(data_Met,'Normalization', 'pdf', 'FaceColor', [0.6350 0.0780 0.1840], 'FaceAlpha', 0.4, 'NumBins', numBins);
[f, xi] = ksdensity(data_N2);
plot(xi, f, 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
[f, xi] = ksdensity(data_Met);
plot(xi, f, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
%ksdensity(data_N2);
%ksdensity(data_Met);
title('Median distances distribution/probability density function ');
legend('Wild Type', '-H3K9me3');
xlabel('Distance');
ylabel('Density');
xlim([0.2 1.2])
ylim([0 5])
hold off;

figure
numBins = 10;
histogram(data_N2, 'Normalization', 'pdf','FaceColor', [0 0.4470 0.7410], 'NumBins', numBins);
hold on;
%histogram(data_Met,'Normalization', 'pdf', 'FaceColor', [0.6350 0.0780 0.1840], 'FaceAlpha', 0.4, 'NumBins', numBins);
histogram(data_Con,'Normalization', 'pdf', 'FaceColor', [0.1660 0.5740 0.1880], 'FaceAlpha', 0.3, 'NumBins', numBins);
%[f, xi] = ksdensity(data_Met);
%plot(xi, f, 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2);
[f, xi] = ksdensity(data_N2);
plot(xi, f, 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
[f, xi] = ksdensity(data_Con);
plot(xi, f, 'Color', [0.1660 0.5740 0.0880], 'LineWidth', 2);
%ksdensity(data_N2);
%ksdensity(data_Met);
title('Median distances distribution/probability density function ');
legend('Wild Type', '(-)H3K9me', '(-)Condensin I');
xlabel('Distance');
ylabel('Density');
xlim([0.2 1.2])
ylim([0 5])
hold off;

%%
TotalTADNum= length(N2DisListAll);

Pval = ones(TotalTADNum, TotalTADNum);

for i = 1:TotalTADNum
    for j = i + 1:TotalTADNum
        % The ttest2 function either computes the two sample t-test or
        % Welch's t-test. The only (real) difference between these two
        % tests is that Welch's t-test does not assume equal variances
        % between the two datasets. Hence we declare 'unequal' in the call.

        % If we set 'Tail' to 'left' we are assuming that N2 < MET2.
        % If we set 'Tail' to 'both' we are assuming that mean(N2
        [~, p, ~, ~] = ttest2(N2DisListAll{i, j}, CONDisListAll{i, j}, 'Tail', 'left', 'Vartype', 'unequal');
        Pval(i, j) = p;
        Pval(j, i) = p;
    end
end
%save('Pval_N2vsMet2.mat','Pval');

%% Create heatmap figure


% close all
figure
imagesc(Pval)
colorbar
axis square
title(['Significant mean changes N2 < Condensin I']);
colormap bone
caxis([0 0.05]);
yticks(1:5:TotalTADNum)
yticklabels({axValues(1:5:TotalTADNum)})
xticks(1:5:TotalTADNum)
xticklabels({axValues(1:5:TotalTADNum)})
xtickangle(90)
%saveas(gcf,'significant_mean_changes_N2_left_MET2','fig');
%print('Anova bin1toLate p05','-dsvg');
%print('Anova bin1toLate p01','-dsvg');
%saveas(figure(6000),'Anova bin1toLate p01','fig');
