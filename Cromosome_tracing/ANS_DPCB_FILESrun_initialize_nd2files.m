% This script opens and converts files from nd2 format to multidemensional matrices for downstream
% processing
% ------------------------------------------------------------------------------------------------------------------------------------------
% Requirements: bfmatlab installed and in the path
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: path to .nd2 file with multiple series(FOV/samples) and illumination channels
%
% series: FOV/position
%
% zstacks for sequential FISH hybs are numbered starting at 3 in increments
% of 1
% hyb1= Channel560,647_seq0003
% hyb2= Channel560,647_seq0004
% hyb3= Channel560,647_seq0005
%
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: Multidemensional matrices for each channel, series, and hyb#
% ------------------------------------------------------------------------------------------------------------------------------------------
% Modified from Ahilya N. Sawh, PhD
% caxis([ 0 1800])for chanel 1
% 17.01.2019
% Version 1.0
%% ------------------------------------------------------------------------------------------------------------------------------------------
clear all
close all
tic

mkdir sequential

bfmatlabpath ='/scicore/home/mangos/pulido0000/Tracing/bfmatlab_5.9.2';
pathString = genpath(bfmatlabpath);
addpath(pathString)

NumSeries = 39;
NumHybs = 28; %the last number of the file
sz = 1608;

%primary probe file

for Hyb = 0
    Seq = Hyb*1;

    data = bfopen(['Channel405,488,560,647_Seq0000Denoised.nd2']);

    for series = 1:NumSeries
        ImageStack = zeros(sz,sz,604);
        for i = 1:604 %planes are organized by z position first then channel
            ImageStack(:,:,i) = data{series, 1}{i,1};
        end

        %channel 1 = 405
        ImageStack405 = zeros(sz,sz,151);
        ImageStack405 = ImageStack(:,:,1:4:604);
        ImageMax = max(ImageStack405,[],3);
%         figure
%         imagesc(ImageMax)
%         axis equal
%         colormap gray
%         title(['405_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        %saveas(gcf,['sequential/405_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        save( ['sequential/405_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack405','-v7.3');
        clear ImageMax ImageStack405

        %channel 2 = 488
        ImageStack488 = zeros(sz,sz,151);
        ImageStack488 = ImageStack(:,:,2:4:604);
        ImageMax = max(ImageStack488,[],3);
%         figure
%         imagesc(ImageMax)
%         axis equal
%         colormap gray
%         title(['488_Hyb' num2str(Hyb) '_FOV' num2str(series)])
%         saveas(gcf,['sequential/488_Hyb' num2str(Hyb) '_FOV' num2str(series)])
         save( ['sequential/488_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack488');
        clear ImageMax ImageStack488

        %channel 3 = 560
        ImageStack560 = zeros(sz,sz,151);
        ImageStack560 = ImageStack(:,:,3:4:604);
        ImageMax = max(ImageStack560,[],3);
%         figure
%         imagesc(ImageMax)
%         axis equal
%         colormap gray
%         title(['560_Hyb' num2str(Hyb) '_FOV' num2str(series)])
%         saveas(gcf,['sequential/560_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        save( ['sequential/560_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack560');
        clear ImageMax ImageStack560

        %channel 4 = 647
        ImageStack647 = zeros(sz,sz,151);
        ImageStack647 = ImageStack(:,:,4:4:604);
        ImageMax = max(ImageStack647,[],3);
%         figure
%         imagesc(ImageMax)
%         axis equal
%         colormap gray
%         title(['647_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        %saveas(gcf,['sequential/647_Hyb' num2str(Hyb) '_FOV' num2str(series)])
        save( ['sequential/647_' num2str(Hyb) '_' num2str(series) '.mat'],'ImageStack647');
        clear ImageMax ImageStack647

    end
    clear data
end


%% secondary 1 to 13 hyb files

for Hyb = 3:15; %3 insted of 1
    Seq = Hyb;
    if Seq <= 9 %change this one for 9 insted of 10
        data = bfopen(['Channel560,488_Seq000' num2str(Seq) '.nd2']);
    elseif Seq >= 10
        data = bfopen(['Channel560,488_Seq00' num2str(Seq) '.nd2']);
    end
    for series = 1:NumSeries
        ImageStack = zeros(sz,sz,302);
        for i = 1:302 %planes are organized by z position first then channel
            ImageStack(:,:,i) = data{series, 1}{i,1};
        end

        %channel 1 = 560
        ImageStack560 = zeros(sz,sz,151);
        ImageStack560 = ImageStack(:,:,1:2:302);
        ImageMax = max(ImageStack560,[],3);
%         figure
%         imagesc(ImageMax)
%         axis equal
%         colormap gray
%         %caxis([ 0 1800])
%         title(['560_Hyb' num2str(Hyb-14) '_FOV' num2str(series)])
        %saveas(gcf,['sequential/560_Hyb' num2str(Hyb-14) '_FOV' num2str(series)])
        save( ['sequential/560_' num2str(Hyb-3) '_' num2str(series) '.mat'],'ImageStack560');
        clear ImageMax ImageStack560

        %channel 2 = 488
        ImageStack488 = zeros(sz,sz,151);
        ImageStack488 = ImageStack(:,:,2:2:302);
        ImageMax = max(ImageStack488,[],3);
%         figure
%         imagesc(ImageMax)
%         axis equal
%         colormap gray
%         title(['488_Hyb' num2str(Hyb-14) '_FOV' num2str(series)])
        %saveas(gcf,['sequential/488_Hyb' num2str(Hyb-14) '_FOV' num2str(series)])
        save( ['sequential/488_' num2str(Hyb-3) '_' num2str(series) '.mat'],'ImageStack488');
        clear ImageMax ImageStack488

        %            %channel 3 = 647
        %             ImageStack647 = zeros(sz,sz,151);
        %             ImageStack647 = ImageStack(:,:,2:2:302);
        %             ImageMax = max(ImageStack647,[],3);
        %             figure
        %             imagesc(ImageMax)
        %             axis equal
        %             colormap gray
        %             title(['647_Hyb' num2str(Hyb+13) '_FOV' num2str(series)])
        %             saveas(gcf,['sequential/647_Hyb' num2str(Hyb+13) '_FOV' num2str(series)])
        %             save( ['sequential/647_' num2str(Hyb+13) '_' num2str(series) '.mat'],'ImageStack647');
        %             clear ImageMax ImageStack647 ImageStack
        %
        close all
    end
    clear data
end

close all

%% secondary 14 to 28 hyb files

for Hyb = 3:17; %3 insted of 1
    Seq = Hyb;
    if Seq <= 9 %change this one for 9 insted of 10
        data = bfopen(['Channel560,488,647_Seq000' num2str(Seq) '.nd2']);
    elseif Seq >= 10
        data = bfopen(['Channel560,488,647_Seq00' num2str(Seq) '.nd2']);
    end
    for series = 1:NumSeries
        ImageStack = zeros(sz,sz,453);
        for i = 1:453 %planes are organized by z position first then channel
            ImageStack(:,:,i) = data{series, 1}{i,1};
        end

        %channel 1 = 560
        ImageStack560 = zeros(sz,sz,151);
        ImageStack560 = ImageStack(:,:,1:3:453);
        ImageMax = max(ImageStack560,[],3);
%         figure
%         imagesc(ImageMax)
%         axis equal
%         colormap gray
%         %caxis([ 0 1800])
%         title(['560_Hyb' num2str(Hyb+14) '_FOV' num2str(series)])
        %saveas(gcf,['sequential/560_Hyb' num2str(Hyb+14) '_FOV' num2str(series)])
        save( ['sequential/560_' num2str(Hyb+11) '_' num2str(series) '.mat'],'ImageStack560');
        clear ImageMax ImageStack560

        %channel 2 = 488
        ImageStack488 = zeros(sz,sz,151);
        ImageStack488 = ImageStack(:,:,2:3:453);
        ImageMax = max(ImageStack488,[],3);
%         figure
%         imagesc(ImageMax)
%         axis equal
%         colormap gray
%         title(['488_Hyb' num2str(Hyb+14) '_FOV' num2str(series)])
%         saveas(gcf,['sequential/488_Hyb' num2str(Hyb+14) '_FOV' num2str(series)])
        save( ['sequential/488_' num2str(Hyb+11) '_' num2str(series) '.mat'],'ImageStack488');
        clear ImageMax ImageStack488

        %channel 3 = 647
        ImageStack647 = zeros(sz,sz,151);
        ImageStack647 = ImageStack(:,:,3:3:453);
        ImageMax = max(ImageStack647,[],3);
%         figure
%         imagesc(ImageMax)
%         axis equal
%         colormap gray
%         title(['647_Hyb' num2str(Hyb+14) '_FOV' num2str(series)])
        %saveas(gcf,['sequential/647_Hyb' num2str(Hyb+14) '_FOV' num2str(series)])
        save( ['sequential/647_' num2str(Hyb+11) '_' num2str(series) '.mat'],'ImageStack647');
        clear ImageMax ImageStack647 ImageStack

        close all
    end
    clear data
end

toc
