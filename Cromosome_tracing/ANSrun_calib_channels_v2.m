% This script calbrates the two channels (560 and 647) used for imaging DNA
% lociin sequential FISH
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: multidimensional matrices (3D image stacks) of Tetraspeck beads 
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: Warping matrix from the 647 channel to the 560 channel (tform + DeltaZ)
% ------------------------------------------------------------------------------------------------------------------------------------------
% Ahilya N. Sawh, PhD
% 23.08.2019
% Version 1.0
% 
% using fitFoci and selectROI functions from the Zhuang lab
%%
clear all
close all

bfmatlabpath ='/scicore/home/mangos/pulido0000/Tracing/bfmatlab_5.9.2';
pathString = genpath(bfmatlabpath);
addpath(pathString)

data = bfopen('220125.nd2');
NumSeries = 1;
sz = 1608;

for series = 1:NumSeries
    ImageStack = zeros(sz,sz,264);
    for i = 1:264 %planes are organized by z position first then channel
        ImageStack(:,:,i) = data{series, 1}{i,1};
    end
    
    
    %channel 1 = 560
    ImageStack560 = zeros(sz,sz,66);
    ImageStack560 = ImageStack(:,:,3:4:264);
    ImageMax = max(ImageStack560,[],3);
    figure
    imagesc(ImageMax)
    axis equal
    colormap gray
    title(['560_calib' num2str(series)])
    saveas(gcf,['560_calib' num2str(series)])
    save( ['560_calib_' num2str(series) '.mat'],'ImageStack560');
    
    
    %channel 2 = 647
    ImageStack647 = zeros(sz,sz,66);
    ImageStack647 = ImageStack(:,:,4:4:264);
    ImageMax = max(ImageStack647,[],3);
    figure
    imagesc(ImageMax)
    axis equal
    colormap gray
    title(['647_calib' num2str(series)])
    saveas(gcf,['647_calib' num2str(series)])
    save( ['647_calib_' num2str(series) '.mat'],'ImageStack647');
    
    
end
clear data


%%
close all
clear all


bfmatlabpath ='/scicore/home/mangos/pulido0000/Tracing/BasicFunctions';
pathString = genpath(bfmatlabpath);
addpath(pathString)

FileName = '560_calib_1.mat';
load(FileName)

for i = 1:size(ImageStack560, 3)
    Std560(i) = std2(ImageStack560(:,:,i)); %calculate the standard deviation of pixel intensity values for the whole slice
end

FileName = '647_calib_1.mat';
load(FileName)

for i = 1:size(ImageStack647, 3)
    Std647(i) = std2(ImageStack647(:,:,i));
end

%

ImageMax560 = max(ImageStack560,[],3);
ImageMax647 = max(ImageStack647,[],3);

figure 
subplot(1,2,1)
imagesc(ImageMax560)
title('561')
colormap gray
axis square
caxis([100 300])
subplot(1,2,2)
imagesc(ImageMax647)
title('647')
colormap gray
axis square
caxis([100 400])

clear RBG
RGB(:,:,1) = ImageMax560/max(max(ImageMax560))*2;
RGB(:,:,2) = ImageMax647/max(max(ImageMax647))*2;
q = size(ImageStack560);
RGB(:,:,3) = zeros(q(1,1),q(1,2));

figure
imagesc(RGB)
axis square
title('beads before transformation')

% align channels in x-y
figure
imagesc(RGB)
axis square
title('select pairs of foci')
WhetherROI = questdlg('Do you want to select ROIs ?');
if(strcmp(WhetherROI, 'Yes'))
    selectROI;
    lmol_match =[];
    rmol_match =[];
    for i = 1:length(roiList)
        [Xfit, Yfit, Zfit] = fitFoci(ImageStack560, roiList(i),2*i-1,1);
        lmol_match = cat(1, lmol_match, [Xfit, Yfit, Zfit]);
        [Xfit, Yfit, Zfit] = fitFoci(ImageStack647, roiList(i),2*i,1);
        rmol_match = cat(1, rmol_match, [Xfit, Yfit, Zfit]);
    end
    tform = cp2tform(rmol_match(:,1:2), lmol_match(:,1:2),'projective'); %right (647) points are moving
    %tform2 = fitgeotrans(rmol_match(:,1:2), lmol_match(:,1:2),'projective'); 
    save('tform.mat','tform');
end

% to use tform
TransImg = imtransform(RGB(:,:,2), tform, 'XData', [1 q(1,1)], 'Ydata', [1 q(1,1)]);
RGB(:,:,2) = TransImg(1:q(1,1), 1:q(1,1));
figure
imagesc(RGB)
axis square
title('beads after transformation')

% uncomment following to use tform2
% TransImg2 = imwarp(RGB(:,:,2), tform2,'OutputView',imref2d(size(RGB)));
% RGB(:,:,2) = TransImg2(1:q(1,1), 1:q(1,1));
% figure
% imagesc(RGB)
% axis square
% title('beads after transformation2')

% align channels in z
GaussEqu = 'a*exp(-(x-b)^2/2/c^2)+d';
%a is the height of the curve's peak, b is the position of the center of the peak 
%and c (the standard deviation, sometimes called the Gaussian RMS width) controls the width of the "bell"

StartPoint = [max(Std560) 9 15 0];
f_left = fit([1:length(Std560)]', Std560', GaussEqu, 'Start', StartPoint);

StartPoint = [max(Std647) 9 15 0];
f_right = fit([1:length(Std647)]', Std647', GaussEqu, 'Start', StartPoint);

figure
subplot(1,2,1)
plot(f_left, 1:length(Std560), Std560);
title('STD 560')
xlabel('Z position')
legend('off');
subplot(1,2,2)
plot(f_right, 1:length(Std647), Std647);
legend('off');
title('STD 647')
xlabel('Z position')

DeltaZ = f_right.b-f_left.b; %position of center peak in two channels 647 minus 560
save('DeltaZ.mat','DeltaZ');




