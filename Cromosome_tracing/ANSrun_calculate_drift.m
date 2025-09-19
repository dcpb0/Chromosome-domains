% This script calculates the sample drift parameters for each sample/FOV in
% sequential FISH by fitting a Gaussian curve to the bead signal in each
% round 
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: multidimensional matrices (3D image stacks) of beads in each round
% of sequential FISH
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: DriftParams file for each sample/FOV containing drift in x,y,z
% from Hyb0
% ------------------------------------------------------------------------------------------------------------------------------------------
% Ahilya N. Sawh, PhD
% 06.03.2019
% Version 1.0
% based on the original script obtained from Siyuan Wang in Zhuang lab run_bead,
% using fitFoci and selectROI functions
%% ------------------------------------------------------------------------------------------------------------------------------------------
tic

clear all
close all

SampleNum = 17;
NumHybs = 28;

FileName = ['sequential/488_0_' num2str(SampleNum) '.mat'];

load(FileName)

ImageMean = max(ImageStack488,[],3);
figure(1)
imagesc(ImageMean)
title(SampleNum)%i added this 
colormap summer %gray i
caxis([100 200])
axis equal

selectROI % choose one bead only

figure(2)
imagesc(ImageMean)
title(SampleNum)%i added this 
colormap gray
axis equal
hold on

for i = 0:NumHybs
    FileName = ['sequential/488_' num2str(i) '_' num2str(SampleNum) '.mat'];
    load(FileName)
    [Xfit(i+1), Yfit(i+1), Zfit(i+1)] = fitFoci(ImageStack488, roiList, i+1,1);
    %saveas(gcf,['DriftParams' num2str(SampleNum)])
    figure(2)
    plot(Xfit(i+1), Yfit(i+1), 'x');
end

hold off

for i = 1:NumHybs
    Xdrift(i) = Xfit(i+1) - Xfit(1);
    Ydrift(i) = Yfit(i+1) - Yfit(1);
    Zdrift(i) = Zfit(i+1) - Zfit(1);
end

%save(['DriftParams' num2str(SampleNum) '.mat'], 'Xdrift', 'Ydrift', 'Zdrift');

toc  


