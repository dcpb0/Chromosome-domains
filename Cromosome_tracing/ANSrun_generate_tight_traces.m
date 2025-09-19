% This script traces chromosomes by nearest neighbour method, using
% watershed segmentation of nuclei and primary probe signal 
% ------------------------------------------------------------------------------------------------------------------------------------------
% Input: multidimensional matrices (3D image stacks) of nuclear signal,
% primary probe signal, DeltaZ, tform, and fitted foci result file for each sample
% 
% create folders before running: traces
% run section by section for each sample,
%
% functions required: traceChromosome_3D_L2, traceChromosome_3D_L1, dis
%
% ------------------------------------------------------------------------------------------------------------------------------------------
% Output: Traces with coordinates of each FOV
% ------------------------------------------------------------------------------------------------------------------------------------------
% modified from Ahilya N. Sawh, PhD by DCPB
% 16.02.2021
% 
% 
%% ------------------------------------------------------------------------------------------------------------------------------------------

tic

clear all
close all

% set parameters
global SampleNum
global flap
global x1
global x2
global y1
global y3

SampleNum = '19'; 
FOV =1; %ROI in the FOV (if number of embryos > 1 in the FOV)

save_traces = 0;

save_age = 1;
man_age = 0;

save_traces
save_age

%for primary probe segmentation

primarytosegment_561 = 0;
primarytosegment_647 = 1;

%%%%%%%%%%maxintensitythreshs%%%%%%%%%%%%%%%%%%%%%%change 405 mainly%%%%%%%

Thresh405 = 90; %
Thresh560 = 170;
Thresh647 = 70; %L233 for the primary probe 
%Thresh488 = 2500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for reconstruction Silvi never changes 
dilateThresh = 5; %
erodeThresh = 5; % L256 

%for binarization and watershed strength%%%change terrmin if traces are too
%many/few per territory

NucbinaryThresh = 70;%L132 %b300 for dapi watersheld, the higher more defined nuclei and the lower bigger 
NucMin= 0.4; %L146

TerrbinaryThresh = 5; %L289 %b15, for primarz, higer less and defiend territores 2:4:30
TerrMin = 0.25; %lower the more territories original 0.25 


TotalTADNum = 43;
TADsToExclude =[]; %DONT FORGET 

Color = colormap(jet(50)); % the number of colors for each foci the number of regions 

%% find and segment nuclei

FileName = ['sequential/405_0_' SampleNum '.mat'];

load(FileName)

ImDAPI = max(ImageStack405,[],3);
figure
imagesc(ImDAPI)
colormap gray
axis equal
caxis([0 1200])
title('Z stack max image')%flattened image of nuclei full xyz stack

% select an embryo FOV in the sample, perform all ds actions in this FOV

selectROI

%% find and segment nuclei II

I = ImageStack405(ceil(y1):ceil(y3),ceil(x1):ceil(x2),:);
figure
ImDAPI = max(I,[],3); 
imagesc(ImDAPI)
colormap gray
axis equal
%caxis([0 900]) %before was 3500
title('selected FOV 405')

%correct for uneven background illumination

se = strel('disk',30); % 60 pixel value, lower removes more background, removes the diference of the value /50
background = imopen(I,se);
figure
imagesc(max(background,[],3))
colormap gray
axis equal
%caxis([0 900]) %before was 3500
title('selected FOV primary probe 405 bkg')

I2 = I - background;
figure
imagesc(max(I2,[],3))
colormap gray
axis equal
%caxis([0 300])
title('selected FOV primary probe 405 minus bkg')

%%  find and segment nuclei III

%set min and max values to even out territory intensity

Thresh405= 150%TO TEST
NucbinaryThresh = 50

I3 = I2;
I3(I3<30) = 0; %80&100 Increase if the DAPI is strong. Inical value 100. %50 
I3(I3>Thresh405) = Thresh405;
figure
imagesc(max(I3,[],3))
colormap gray
axis equal
%caxis([0 500])
title('selected FOV 405 thresh')

%3D binary

 %TO TEST  

BW = imbinarize(I3, NucbinaryThresh);
nucBW = imfill(BW, 'holes');
figure
imshowpair(max(I3,[],3),max(nucBW,[],3),'montage')
%title('before and after 3D binary')
title(['before and after 3D binary', ' NucbiTh ' , num2str(NucbinaryThresh)])

%3D dist transform
D = -bwdist(~nucBW); %This is because the image is binary  computes the Euclidean distance transform of the binary image BW. For each pixel in BW, the distance transform assigns a number that is the distance between that pixel and the nearest nonzero pixel of BW.

%set background to its own catchment basin - i.e. where BW is black
D(~nucBW) = -Inf;

%3D watershed
%NucMin= 0.4;%TO TEST 
D = imhmin(D,NucMin); %the height threshold for suppressing shallow minima
L = watershed(D);
figure
imshow(label2rgb(L(:,:,50),'jet','w'))
title(['watershed NucMin', ' NucbiTh ' , num2str(NucbinaryThresh)])
%title('watershed NucMin')%Added

L = imdilate(L,true(5));
nucL = L-1;
clear L

Lmax = max(nucL,[],3);
nucLrgb = label2rgb(Lmax,'jet','w','shuffle');
figure
imshow(nucLrgb)
%title('colored watershed label matrix nuclei')
title(['colored watershed label matrix nuclei', ' NucbiTh ' , num2str(NucbinaryThresh)])


figure
imagesc(ImDAPI)
colormap gray
hold on
himage = imshow(nucLrgb);
himage.AlphaData = 0.2;
%title('nuclei labels superimposed transparently on DAPI')
title(['nuclei labels superimposed transparently on DAPI', ' NucbiTh ' , num2str(NucbinaryThresh)])


age = max(max(max(nucL)))

nstats = regionprops3(nucL, 'all');ncenters = nstats.Centroid;

%3D nuclear volume

% figure 
% isosurface(nucL, 0.5) 
% axis equal
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('3D territory volumes with all foci in FOV zeroed')


%% load drift calibration 

load('DeltaZ.mat');
load('tform.mat');

%% segment based on primary probe 647

if primarytosegment_647 == 1
    FileName = ['sequential/647_0_' SampleNum];
    load(FileName)
    TransformOrNot = 1;
    
if TransformOrNot == 1
    for m = 1:size(ImageStack647,3)
        ImageStack647(:,:,m) = imtransform(ImageStack647(:,:,m), tform, 'XData', [1 1608], 'Ydata', [1 1608]);
    end
end

Im647 = max(ImageStack647,[],3);
figure
imagesc(Im647)
colormap gray
axis equal
title('z stack max image 647 primary ')
caxis([100 400])

%% segment based on primary probe 647 II

Thresh647 =200;%TO TEST

I = ImageStack647(ceil(y1):ceil(y3),ceil(x1):ceil(x2),:);

%correct for uneven background illumination

Im647 = max(I,[],3);
figure
imagesc(Im647)
colormap gray
axis equal
caxis([0 1800])%B500
title('selected FOV primary probe 647')

se = strel('disk',10); %10 is the vale of the pixel, denoised
background = imopen(I,se);
figure
imagesc(max(background,[],3))
colormap gray
axis equal
caxis([0 800])%B500
title('selected FOV primary probe 647 bkg')

I2 = I - background;
figure
imagesc(max(I2,[],3))
colormap gray
axis equal
caxis([0 300])
title('selected FOV primary probe 647 minus bkg')

%%set min and max values to even out territory intensity


I3 = I2;
I3(I3<10) = 0;%10 lower when primary is weak 
I3(I3>Thresh647) = Thresh647;
figure
imagesc(max(I3,[],3))
colormap gray
axis equal
caxis([0 100])
title('selected FOV primary probe 647 thresh')


% refine by dilate+erode
% se = strel('disk', dilateThresh);
% Id = imdilate(I3, se);
% Icbr = imreconstruct(imcomplement(Id), imcomplement(I3));
% Icbr = imcomplement(Icbr);
% figure
% imagesc(max(Icbr,[],3))
% colormap gray
% axis equal
% title('closing by reconstruction');

%erodeThresh = 5 %added to test in the pictures
se = strel('disk', erodeThresh);
%Icbre = imerode(Icbr, se);
Icbre = imerode(I3, se);
%Icbrobr = imreconstruct(Icbre, Icbr);
Icbrobr = imreconstruct(Icbre, I3);
figure
imagesc(max(Icbrobr,[],3))
colormap gray
axis equal
title('closing and opening by reconstruction');

% remove noise
% J = medfilt3(Icbrobr, [5 5 5]);
% figure
% imshowpair(max(Icbrobr,[],3),max(J,[],3),'montage')
% title('before and after median filter')
% Icbrobr = J;

% binarize 
%2D binary
% BW = imbinarize(Icbrobr(:,:,45), TerrbinaryThresh);
% BW = imfill(BW, 'holes');
% radius = 10;
% decomposition = 0;
% se = strel('disk', radius, decomposition);
% BW = imclose(BW, se);
% figure
% imshowpair(Icbrobr(:,:,45),BW,'montage')
% title('before and after 2Dbinary')

%% segment based on primary probe 647 III

%3D binary
TerrbinaryThresh = 60 %TO TEST
TerrMin = 0.2

BW = imbinarize(Icbrobr, TerrbinaryThresh);
BW = imfill(BW, 'holes');
BW = imclose(BW, se);
figure
imshowpair(max(Icbrobr,[],3),max(BW,[],3),'montage')
%title(['before and after 3Dbinary', ' TerrBiTh ' , num2str(TerrbinaryThresh)])
%title('before and after 3Dbinary')

%3D dist transform
D = -bwdist(~BW); 

%set background to its own catchment basin - where BW is black
D(~BW) = -Inf;

%3D watershed
%TerrMin = 0.25 %TO TESt 0.25
D = imhmin(D,TerrMin); %the height threshold for suppressing shallow minima
L = watershed(D);
% figure
% imshow(label2rgb(L(:,:,50),'jet','w'))

L = imdilate(L,true(5));
L = L-1; %remove background catchment basin - set to 0

L(find(~nucL))=0;

Lmax = max(L,[],3);
Lrgb = label2rgb(Lmax,'jet','w','shuffle');
figure
imshow(Lrgb)
title(['colored watershed label matrix territory 647', ' TerrBiTh ' , num2str(TerrbinaryThresh)])
%title('colored watershed label matrix territory 647')

figure
imagesc(Im647)
colormap gray
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.2;
himage = imshow(nucLrgb);
himage.AlphaData = 0.05;
title(['territory labels superimposed transparently on 647', ' TerrBiTh ' , num2str(TerrbinaryThresh)])
%title('territory labels superimposed transparently on 647')

numterr = max(max(max(L)))
% 
% figure
% imagesc(ImDAPI)
% colormap gray
% hold on
% himage = imshow(Lrgb);
% himage.AlphaData = 0.4;
% himage = imshow(nucLrgb);
% himage.AlphaData = 0.07;
% %title('nuclei and territory labels superimposed transparently on DAPI')
% title(['nuclei and territory labels superimposed transparently on DAPI', ' TerrBiTh ' , num2str(TerrbinaryThresh)])

L2 = L;
L2max = max(L2,[],3);
L2rgb = label2rgb(L2max,'jet','w','shuffle');
% figure
% imshow(L2rgb)
% title('cells territories 647')
% axis equal

% clear L

% % here split nucL into pcell or non-pcell 
% 
% pnucL = zeros(size(nucL));
% pnucL(find(nucL == pcellnucID)) = 1; %copy coords of pcellnucleus
% pnucLmax = max(pnucL,[],3);
% % figure
% % imagesc(pnucLmax)
% % axis equal
% % title('p nucleus')
% 
% somnucL = nucL;
% somnucL(find(nucL == pcellnucID)) = 0; %remove the pcellnucleus
% somnucLmax = max(somnucL,[],3);
% % figure
% % imagesc(somnucLmax)
% % axis equal
% % title('somatic nuclei')

% %%here split L2 into pcell or non-pcell 
% 
% pL2 = L2;
% pL2(find(somnucL))= 0;  
% pL2max = max(pL2,[],3);
% pL2rgb = label2rgb(pL2max,'jet','w','shuffle');
% figure
% imshow(pL2rgb)
% title('pcell territories 647')
% axis equal
% 

%modified this part 
% L2 = L2; %l2 removed for l 
% L2(find(pnucL))= 0;  
% L2max = max(L2,[],3);
% L2rgb = label2rgb(L2max,'jet','w','shuffle');
% figure
% imshow(L2rgb)
% title('cells territories 647')
% axis equal

end

%% load and crop results file to FOV

load(['Results/result' num2str(SampleNum) '.mat']);

% warp the 647 image into the 561 channel in z.
for i = 14:2:TotalTADNum-1
   Zfit{i} = Zfit{i}-DeltaZ;
end

if ~isempty(TADsToExclude)
    for i = TADsToExclude 
        Xfit{1,i} = [];
        Yfit{1,i} = [];
        Zfit{1,i} = [];
        Xgof{1,i} = [];
        Ygof{1,i} = [];
        Zgof{1,i} = [];
        Intensity{1,i} = [];
    end
end

%convert to pixel
for i = 1:length(Intensity)
    Yfit{i} = 1608*110/1000-Yfit{i}; 
    Yfit{i} = Yfit{i}/110*1000; %pxl
    Xfit{i} = Xfit{i}/110*1000; %pxl
    Zfit{i} = Zfit{i}/0.2; % z stack no.
end
figure
hold on
for i = 1:length(Intensity)
    scatter3(Xfit{i}, Yfit{i}, Zfit{i}, 'ok', 'MarkerFaceColor', Color (i,:));
    text(Xfit{i}+0.05, Yfit{i}, Zfit{i}, num2str(i), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',12); %write numbers on the image
end
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'YDir','reverse')
title('All foci - pxl')

Xfitnew = cell(1,TotalTADNum);
Yfitnew = cell(1,TotalTADNum);
Zfitnew = cell(1,TotalTADNum);
Xgofnew = cell(1,TotalTADNum);
Ygofnew = cell(1,TotalTADNum);
Zgofnew = cell(1,TotalTADNum);
Intensitynew = cell(1,TotalTADNum);
for i=1:TotalTADNum
    q = 0;
    if ~isempty(Xfit{1,i})
        for j=1:length(Xfit{1,i})
            if Xfit{1,i}(1,j) >= x1 && Xfit{1,i}(1,j) <= x2 && Yfit{1,i}(1,j) >= y1 && Yfit{1,i}(1,j) <= y3 %inside FOV coordinates
                q = q+1;
                Xfitnew{1,i}(1,q) = Xfit{1,i}(1,j);
                Yfitnew{1,i}(1,q) = Yfit{1,i}(1,j);
                Zfitnew{1,i}(1,q) = Zfit{1,i}(1,j);
                Xgofnew{1,i}(1,q) = Xgof{1,i}(1,j);
                Ygofnew{1,i}(1,q) = Ygof{1,i}(1,j);
                Zgofnew{1,i}(1,q) = Zgof{1,i}(1,j);
                Intensitynew{1,i}(1,q) = Intensity{1,i}(1,j);
            end
        end
    end
end

clear Xfit Yfit Zfit Xgof Ygof Zgof Intensity

Xfit = Xfitnew;
Yfit = Yfitnew;
Zfit = Zfitnew;
Xgof = Xgofnew;
Ygof = Ygofnew;
Zgof = Zgofnew;
Intensity = Intensitynew;

figure
hold on
for i = 1:length(Intensity)
    scatter3(Xfit{i}, Yfit{i}, Zfit{i}, 'ok', 'MarkerFaceColor', Color (i,:));
    text(Xfit{i}+0.05, Yfit{i}, Zfit{i}, num2str(i), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',12); %write numbers on the image
end
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'YDir','reverse')
title('All foci in FOV')

%zero the foci locations to the cropped FOV in x and y so that they align
%with the cropped image from prev steps

Xfit_new = cell(1,TotalTADNum);
Yfit_new = cell(1,TotalTADNum);
for i=1:TotalTADNum
    for j=1:length(Xfit{1,i})
        Xfit_new{1,i} = Xfit{1,i} - x1;
        Yfit_new{1,i} = Yfit{1,i} - y1;
    end
end

figure
imagesc(Im647)%we added this with ANS to see the FOCI in the primary figure 
colormap gray
hold on
for i = 1:length(Intensity)
    scatter3(Xfit_new{i}, Yfit_new{i}, Zfit{i}, 'ok', 'MarkerFaceColor', Color (i,:));
    text(Xfit{i}+0.05, Yfit{i}, Zfit{i}, num2str(i), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',12); %write numbers on the image
end
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'YDir','reverse')
title('All foci in FOV zeroed')

% %This makes a figure with the volume of the terriroty in 3D with all the foci
% 
figure 
isosurface(L2, 0.5) 
hold on
for i = 1:length(Intensity)
    scatter3(Xfit_new{i}, Yfit_new{i}, Zfit{i}, 'ok', 'MarkerFaceColor', Color (i,:));
    text(Xfit{i}+0.05, Yfit{i}, Zfit{i}, num2str(i), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',12); %write numbers on the image
end
hold off
axis equal
xlabel('x');
ylabel('y');
zlabel('z');
title('3D territory volumes with all foci in FOV zeroed')


%% trace chromosomes inside 647 territories 

clear pterrIDs
clear somterrIDs

somterrIDs = unique(L2(:));
idx = unique(L2(:)) >0;
somterrIDs = somterrIDs(idx);

figure
imagesc(Im647)
colormap gray
hold on
himage = imshow(L2rgb);
himage.AlphaData = 0.25;
himage = imshow(nucLrgb);
himage.AlphaData = 0.1;

n = 0;
for i = 1:length(somterrIDs)
    n = n+1;
    w = somterrIDs(n,1);
    [Trace, TAD_id] = traceChromosome_3D_L2(L2, w, Xfit_new, Yfit_new, Zfit, Intensity);
    if flap ~= 1
        if numel(Trace) == 1
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 2
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 3
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'dg', 'MarkerFaceColor', 'g')
            plot3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'-g');
            text(Trace{1,3}(:,1), Trace{1,3}(:,2),Trace{1,3}(:,3), num2str(TAD_id{3}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 4
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'dg', 'MarkerFaceColor', 'g')
            plot3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'-g');
            text(Trace{1,3}(:,1), Trace{1,3}(:,2),Trace{1,3}(:,3), num2str(TAD_id{3}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,4}(:,1), Trace{1,4}(:,2), Trace{1,4}(:,3),'^k', 'MarkerFaceColor', 'k')
            plot3(Trace{1,4}(:,1), Trace{1,4}(:,2), Trace{1,4}(:,3),'-k');
            text(Trace{1,4}(:,1), Trace{1,4}(:,2),Trace{1,4}(:,3), num2str(TAD_id{4}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        end
    end
    if save_traces == 1 && save_age ==0
        save(['traces/' num2str(man_age) '_cell_Traces_Sample' SampleNum '_FOV' num2str(FOV) '_roi647_' num2str(w) '.mat'],...
            'Trace', 'TAD_id', 'L2', 'nucL', 'TADsToExclude', 'x1', 'x2', 'y1', 'y3');
    end
    if save_traces == 1 && save_age ==1
        save(['traces/' num2str(age) '_cell_Traces_Sample' SampleNum '_FOV' num2str(FOV) '_roi647_' num2str(w) '.mat'],...
            'Trace', 'TAD_id', 'L2', 'nucL', 'TADsToExclude', 'x1', 'x2', 'y1', 'y3');
    end
end
title('traces on segmentation stain 647')
axis equal
hold off

%% segment based on primary probe 561

if primarytosegment_561 == 1
    FileName = ['sequential/560_0_' SampleNum];
    load(FileName)
    TransformOrNot = 0; 


    if TransformOrNot == 1
        for m = 1:size(ImageStack560,3)
            ImageStack560(:,:,m) = imtransform(ImageStack560(:,:,m), tform, 'XData', [1 1608], 'Ydata', [1 1608]);
        end
    end

    Im560 = max(ImageStack560,[],3);
    figure
    imagesc(Im560)
    colormap gray
    axis equal
    title('z stack max image 561 primary ')
    caxis([100 700])

    I = ImageStack560(ceil(y1):ceil(y3),ceil(x1):ceil(x2),:);

    %correct for uneven background illumination
    Im560 = max(I,[],3);
    figure
    imagesc(Im560)
    colormap gray
    axis equal
    caxis([0 1000])
    title('selected FOV primary probe 560')

    se = strel('disk',10);
    background = imopen(I,se);
    figure
    imagesc(max(background,[],3))
    colormap gray
    axis equal
    caxis([0 1000])
    title('selected FOV primary probe 560 bkg')

    I2 = I - background;
    figure
    imagesc(max(I2,[],3))
    colormap gray
    axis equal
    caxis([0 500])
    title('selected FOV primary probe 560 minus bkg')

    %set min and max values to even out territory intensity

    I3 = I2;
    I3(I3<50) = 0;
    I3(I3>Thresh560) = Thresh560;
    figure
    imagesc(max(I3,[],3))
    colormap gray
    axis equal
    caxis([0 400])
    title('selected FOV primary probe 530 thresh')

    % refine by dilate+erode
    % se = strel('disk', dilateThresh);
    % Id = imdilate(I3, se);
    % Icbr = imreconstruct(imcomplement(Id), imcomplement(I));
    % Icbr = imcomplement(Icbr);
    % figure
    % imagesc(max(Icbr,[],3))
    % colormap gray
    % axis equal
    % title('closing by reconstruction');

    se = strel('disk', erodeThresh);
    %Icbre = imerode(Icbr, se);
    Icbre = imerode(I3, se);
    %Icbrobr = imreconstruct(Icbre, Icbr);
    Icbrobr = imreconstruct(Icbre, I3);
    figure
    imagesc(max(Icbrobr,[],3))
    colormap gray
    axis equal
    title('closing and opening by reconstruction');

    % remove noise
    % J = medfilt3(Icbrobr, [5 5 5]);
    % figure
    % imshowpair(max(Icbrobr,[],3),max(J,[],3),'montage')
    % title('before and after median filter')
    % Icbrobr = J;

    % binarize 
    %2D binary
    % BW = imbinarize(Icbrobr(:,:,45), TerrbinaryThresh);
    % BW = imfill(BW, 'holes');
    % radius = 10;
    % decomposition = 0;
    % se = strel('disk', radius, decomposition);
    % BW = imclose(BW, se);
    % figure
    % imshowpair(Icbrobr(:,:,45),BW,'montage')
    % title('before and after 2Dbinary')

    %3D binary
    BW = imbinarize(Icbrobr, TerrbinaryThresh);
    BW = imfill(BW, 'holes');
    BW = imclose(BW, se);
    figure
    imshowpair(max(Icbrobr,[],3),max(BW,[],3),'montage')
    title('before and after 3Dbinary')

    %3D dist transform
    D = -bwdist(~BW); 

    %set background to its own catchment basin - where BW is 0
    D(~BW) = -Inf;

    %3D watershed
    D = imhmin(D,TerrMin); %the height threshold for suppressing shallow minima
    L = watershed(D);
    figure
    imshow(label2rgb(L(:,:,50),'jet','w'))
    title('watershed slice 50')

    % figure
    % heatmap(L(:,:,50))
    % title('watershed slice 50')

    L = imdilate(L,true(5));
    L = L-1; %remove background catchment basin - set to 0

    %%if pixels were bkg in nuclei label matrix and
    %%if yes convert those pixels to bkg in territory label matrix (zero)

    L(find(~nucL))=0;

    Lmax = max(L,[],3);
    Lrgb = label2rgb(Lmax,'jet','w','shuffle');
    imshow(Lrgb)
    title('colored watershed label matrix territory 560')

    figure
    imagesc(Im560)
    colormap gray
    hold on
    himage = imshow(Lrgb);
    himage.AlphaData = 0.3;
    title('territory labels superimposed transparently on 560')

    numterr = max(max(max(L)))

    L1 = L;
    clear L 

    % here split nucL into pcell or non-pcell 

    pnucL = zeros(size(nucL));
    pnucL(find(nucL == pcellnucID)) = 1; %copy coords of pcellnucleus
    pnucLmax = max(pnucL,[],3);
    figure
    imagesc(pnucLmax)
    axis equal
    title('p nucleus')

    somnucL = nucL;
    somnucL(find(nucL == pcellnucID)) = 0; %remove the pcellnucleus
    somnucLmax = max(somnucL,[],3);
    figure
    imagesc(somnucLmax)
    axis equal
    title('somatic nuclei')

    %%here split L1 into pcell or non-pcell 

    pL1 = L1;
    pL1(find(somnucL))= 0;  
    pL1max = max(pL1,[],3);
    pL1rgb = label2rgb(pL1max,'jet','w','shuffle');
    figure
    imshow(pL1rgb)
    title('pcell territories 560')
    axis equal

    npL1 = L1;
    npL1(find(pnucL))= 0;  
    npL1max = max(npL1,[],3);
    npL1rgb = label2rgb(npL1max,'jet','w','shuffle');
    figure
    imshow(npL1rgb)
    title('nonpcell territories 560')
    axis equal

end 

%% load and crop results file to FOV


%% trace chromosomes inside 560 territories (L1) split by p and nonp cell

if primarytosegment_561 == 1 %I ADDED THIS DCPB
    
pterrIDs = unique(pL1(:));
idx = unique(pL1(:)) >0;
pterrIDs = pterrIDs(idx);

figure
imagesc(Im560)
colormap gray
hold on
himage = imshow(pL1rgb);
himage.AlphaData = 0.25;
himage = imshow(nucLrgb);
himage.AlphaData = 0.1;

n = 0;
for i = 1:length(pterrIDs)
    n = n+1;
    w = pterrIDs(n,1);
    [Trace, TAD_id] = traceChromosome_3D_L1(pL1, w, Xfit_new, Yfit_new, Zfit, Intensity);
    if flap ~= 1
        if numel(Trace) == 1
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 2
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 3
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'dg', 'MarkerFaceColor', 'g')
            plot3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'-g');
            text(Trace{1,3}(:,1), Trace{1,3}(:,2),Trace{1,3}(:,3), num2str(TAD_id{3}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 4
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'dg', 'MarkerFaceColor', 'g')
            plot3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'-g');
            text(Trace{1,3}(:,1), Trace{1,3}(:,2),Trace{1,3}(:,3), num2str(TAD_id{3}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,4}(:,1), Trace{1,4}(:,2), Trace{1,4}(:,3),'^k', 'MarkerFaceColor', 'k')
            plot3(Trace{1,4}(:,1), Trace{1,4}(:,2), Trace{1,4}(:,3),'-k');
            text(Trace{1,4}(:,1), Trace{1,4}(:,2),Trace{1,4}(:,3), num2str(TAD_id{4}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        end
    end
    if save_traces == 1 && save_age ==0
        save(['traces_p\' num2str(man_age) '_cell_Traces_Sample' SampleNum '_FOV' num2str(FOV) '_roi561_' num2str(w) '.mat'],...
            'Trace', 'TAD_id', 'pL1', 'nucL', 'TADsToExclude', 'x1', 'x2', 'y1', 'y3');
    end
    if save_traces == 1 && save_age ==1
        save(['traces_p\' num2str(age) '_cell_Traces_Sample' SampleNum '_FOV' num2str(FOV) '_roi561_' num2str(w) '.mat'],...
            'Trace', 'TAD_id', 'pL1', 'nucL', 'TADsToExclude', 'x1', 'x2', 'y1', 'y3');
    end
end
title('P traces on segmentation stain 561')
axis equal
hold off


somterrIDs = unique(npL1(:));
idx = unique(npL1(:)) >0;
somterrIDs = somterrIDs(idx);

figure
imagesc(Im560)
colormap gray
hold on
himage = imshow(npL1rgb);
himage.AlphaData = 0.25;
himage = imshow(nucLrgb);
himage.AlphaData = 0.1;

n = 0;
for i = 1:length(somterrIDs)
    n = n+1;
    w = somterrIDs(n,1);
    [Trace, TAD_id] = traceChromosome_3D_L1(npL1, w, Xfit_new, Yfit_new, Zfit, Intensity);
    if flap ~= 1
        if numel(Trace) == 1
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 2
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 3
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'dg', 'MarkerFaceColor', 'g')
            plot3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'-g');
            text(Trace{1,3}(:,1), Trace{1,3}(:,2),Trace{1,3}(:,3), num2str(TAD_id{3}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        elseif numel(Trace) == 4
            scatter3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'or', 'MarkerFaceColor', 'r')
            plot3(Trace{1,1}(:,1), Trace{1,1}(:,2), Trace{1,1}(:,3),'-r');
            text(Trace{1,1}(:,1), Trace{1,1}(:,2),Trace{1,1}(:,3), num2str(TAD_id{1}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'sb', 'MarkerFaceColor', 'b')
            plot3(Trace{1,2}(:,1), Trace{1,2}(:,2), Trace{1,2}(:,3),'-b');
            text(Trace{1,2}(:,1), Trace{1,2}(:,2),Trace{1,2}(:,3), num2str(TAD_id{2}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'dg', 'MarkerFaceColor', 'g')
            plot3(Trace{1,3}(:,1), Trace{1,3}(:,2), Trace{1,3}(:,3),'-g');
            text(Trace{1,3}(:,1), Trace{1,3}(:,2),Trace{1,3}(:,3), num2str(TAD_id{3}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
            scatter3(Trace{1,4}(:,1), Trace{1,4}(:,2), Trace{1,4}(:,3),'^k', 'MarkerFaceColor', 'k')
            plot3(Trace{1,4}(:,1), Trace{1,4}(:,2), Trace{1,4}(:,3),'-k');
            text(Trace{1,4}(:,1), Trace{1,4}(:,2),Trace{1,4}(:,3), num2str(TAD_id{4}), 'Color', 'Black', ...
                'FontWeight', 'Bold', 'FontSize',8);
        end
    end
    if save_traces == 1 && save_age ==0
        save(['traces\' num2str(man_age) '_cell_Traces_Sample' SampleNum '_FOV' num2str(FOV) '_roi561_' num2str(w) '.mat'],...
            'Trace', 'TAD_id', 'npL1', 'nucL', 'TADsToExclude', 'x1', 'x2', 'y1', 'y3');
    end
    if save_traces == 1 && save_age ==1
        save(['traces\' num2str(age) '_cell_Traces_Sample' SampleNum '_FOV' num2str(FOV) '_roi561_' num2str(w) '.mat'],...
            'Trace', 'TAD_id', 'npL1', 'nucL', 'TADsToExclude', 'x1', 'x2', 'y1', 'y3');
    end
end
title('nonP traces on segmentation stain 561')
axis equal
hold off

end 

toc

           