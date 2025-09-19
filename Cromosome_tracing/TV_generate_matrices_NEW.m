%% This script generates the Chr variable and corresponding L2, NucL,Mean,
% Median, STD, SEM, NofData, and DisList matrices.
%
% This script is adapted from Ahylia Sawh's 'ANSrun_generate_matrices.m'
% script. The old script could produce the variables with no issues,
% however speed was a huge drawback. It was reported that running the
% script could take days, or sometimes weeks if trying to load several
% folders at once.
%
% Now we are able to produce the same information in under 10 minutes
% per individual folder (run times can be longer if saving L2 and NucL
% data as well)
%
% TV 06.2023

% clear workspace
clear
close all

% These values can be adjusted depending on individual cases
TotalTADNum = 43;
sz = 1608;
AgeStart = 2; % change these to specify embryo age
AgeStop = 140;
AgeRange = ([ num2str(AgeStart) 'to' num2str(AgeStop) 'cell']);
dirName = AgeRange;

% Check if the user wants to save all of the information, or just to
% generate it into the workspace for storage later.
save_flag = check_save_input();

% Check if the folder exists, if not generate it, if it does, confirm user
% wants to continue
if save_flag == 1
    if ~isfolder(dirName)
        mkdir(dirName)
    else
        check_folder_input();
    end
end

% Check if user wishes to save the Chr, L2, and NucL variables
c_flag = check_C_input();
l_flag = check_L_input();
n_flag = check_N_input();

% start timer
tic;

% Define all folder names and prereqs to read in the files in the Trace
% folders
% DON'T FORGET to initialize: foldernames, files_in_dir, and num_files or
% you will be missing information at the end of compilation
% Include the bracing, even for only one trace
foldernames{1} = 'traces_10';
files_in_dir{1} = dir([foldernames{1} '/*.mat']);
num_files(1) = length(files_in_dir{1});

foldernames{2} = 'traces_15';
files_in_dir{2} = dir([foldernames{2} '/*.mat']);
num_files(2) = length(files_in_dir{2});

% foldernames{3} = 'traces_15';
% files_in_dir{3} = dir([foldernames{3} '/*.mat']);
% num_files(3) = length(files_in_dir{3});

% Run through the folders and check the file names for the correct age
% range, if they are not in the right age range, remove them from the
% variable to decrease the number of files we're processing.
for f_ind = 1:length(foldernames)
    for i = num_files(f_ind):-1:1
        % Extract all strings of numbers within the name of the file at 
        % index i
        [ind1, ind2] = regexp(files_in_dir{f_ind}(i).name, '\d+');
        % The age of our trace is stored in the first element of our index
        % variables
        age = str2double(extractBetween(files_in_dir{f_ind}(i).name, ind1(1), ind2(1)));
        if age < AgeStart || age > AgeStop
            files_in_dir{f_ind}(i) = [];
        end
    end

    % update the number of files to the condensed list of suitable files
    num_files(f_ind) = length(files_in_dir{f_ind}); %#ok<*SAGROW>
end
% Extract the Segment Channel from the filename, all files must have the
% same channel.
SegChannel = str2double(extractBetween(files_in_dir{1}(1).name, ...
    ind1(4), ind2(4)));

% sum total of all relevant traces in all folders
total_num_files = sum(num_files);

% Preallocate storage space for variables (not doing this can cause data
% corruption and severe decrease in performance)
% This is naively assuming that there will be some invalid traces and not
% too many traces with multiple counts
for i = total_num_files:-1:1
    if l_flag == 1
        L2s(i).L2 = -1;
        L2s(i).TraceNo = -1;
        L2s(i).TraceName = ' ';
    end
    if n_flag == 1
        NucLs(i).nucL = -1;
        NucLs(i).TraceNo = -1;
        NucLs(i).TraceName = ' ';
    end
    if c_flag == 1
        Chr(i).x = zeros(TotalTADNum, 1);
        Chr(i).y = zeros(TotalTADNum, 1);
        Chr(i).z = zeros(TotalTADNum, 1);
        Chr(i).r = zeros(TotalTADNum, 1);
        Chr(i).RoG = -10;
        Chr(i).Age = -1;
        Chr(i).TraceNo = -1;
        Chr(i).TracesInTerritory = -1;
        Chr(i).TraceName = ' ';
    end
end

% initialize a progress bar to see progress of reading files
wait_bar = waitbar(0, '1', 'Name', 'Reading Files', ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');

setappdata(wait_bar, 'canceling', 0);

% start a counter (takes place of 'n' variable in old script)
varcount = 0;
numsteps = 0;
% Iterate through all files in our folder
for f_ind = 1:length(foldernames)
    for file_ind = 1:num_files(f_ind)
        % If cancel button on progress bar is pressed, program stops executing
        if getappdata(wait_bar, 'canceling')
            sprintf('User cancelled...exiting...')
            delete(wait_bar)
            return
        end
        % Update progress bar
        prg = numsteps/total_num_files;
        waitbar(prg, wait_bar, sprintf('percentage: %2.4f %%', 100 * prg))
        numsteps = numsteps + 1;

        % Create tracename
        TraceName = [foldernames{f_ind} '/' files_in_dir{f_ind}(file_ind).name];
        % Extract all strings of numbers within the given filename
        [ind1, ind2] = regexp(files_in_dir{f_ind}(file_ind).name, '\d+');
        % The age of our trace is stored in the first element of our index
        % variables
        age = str2double(extractBetween(files_in_dir{f_ind}(file_ind).name, ...
            ind1(1), ind2(1)));

        % load our data
        load(TraceName)

        % Run through every element in our trace (also works if only 1 var)
        for trace_count = 1:length(Trace)
            % Ensure trace is what we're interested in.
            if ~check_prereqs(Trace, TAD_id, trace_count)
                % skip trace if it has multiple foci per TAD/region (remove
                % ambiguous traces) or if it doesn't have information
                continue;
            end
            % Update counter
            varcount = varcount + 1;
            Chr(varcount).x = zeros(TotalTADNum, 1);
            Chr(varcount).y = zeros(TotalTADNum, 1);
            Chr(varcount).z = zeros(TotalTADNum, 1);
            Chr(varcount).r = zeros(TotalTADNum, 1);

            if c_flag == 1
                % Iterate through the trace contents
                for k = 1:length(TAD_id{1, trace_count})
                    % rescale Trace values to convert to um from pixel units
                    Trace{1, trace_count}(k, 1) = Trace{1, trace_count}...
                        (k, 1) * 0.11;
    
                    Trace{1, trace_count}(k, 2) = Trace{1, trace_count}...
                        (k, 2) * 0.11;
    
                    Trace{1, trace_count}(k, 3) = Trace{1, trace_count}...
                        (k, 3) * 0.2;
    
                    % Store x, y, z values in Chr
                    Chr(varcount).x(TAD_id{1, trace_count}(k, 1)) = ...
                        Trace{1, trace_count}(k, 1);
    
                    Chr(varcount).y(TAD_id{1, trace_count}(k, 1)) = sz ...
                        * 0.11 - Trace{1, trace_count}(k, 2);
    
                    Chr(varcount).z(TAD_id{1, trace_count}(k, 1)) = ...
                        Trace{1, trace_count}(k, 3);
    
                    Chr(varcount).r(TAD_id{1, trace_count}(k, 1)) = 1;
                end
                % Store rest of information
                Chr(varcount).RoG = rog(Trace{1, trace_count});
                Chr(varcount).Age = age;
                Chr(varcount).TraceNo = trace_count;
                Chr(varcount).TracesInTerritory = length(Trace);
                Chr(varcount).TraceName = TraceName;
            end

            % Check if we want to continue with L and NucL data
            if l_flag == 1
                if SegChannel == 647
                    L2s(varcount).L2 = L2; % 647 segmentation = L2
                    L2s(varcount).TraceNo = trace_count;
                    L2s(varcount).TraceName = TraceName;
                else
                    L1s(varcount).L1 = L1; % 561 segmentation = L1
                    L1s(varcount).TraceNo = trace_count;
                    L1s(varcount).TraceName = TraceName;
                end
            end
            if n_flag == 1
                NucLs(varcount).nucL = nucL;
                NucLs(varcount).TraceNo = trace_count;
                NucLs(varcount).TraceName = TraceName;
            end
        end
    end
end
% Need to call this otherwise you literally cannot close the progress bar
delete(wait_bar)

% Check for and remove unupdated information
if c_flag == 1
    Chr([Chr.Age] == -1) = [];
end
if l_flag == 1
    L2s([L2s.TraceNo] == -1) = [];
end
if n_flag == 1
    NucLs([NucLs.TraceNo] == -1) = [];
end
toc
%% Calculate Matrices
tic
% Initialize new progress bar
wait_bar = waitbar(0, '1', 'Name', 'Writing matrices', ...
    'CreateCancelBtn', 'setappdata(gcbf,''canceling'',1)');

setappdata(wait_bar, 'canceling', 0);

if c_flag == 1
    % Pre-allocate storage space for the matrices
    Mean = zeros(TotalTADNum, TotalTADNum);
    Std = zeros(TotalTADNum, TotalTADNum);
    SEM = zeros(TotalTADNum, TotalTADNum);
    NofData = zeros(TotalTADNum, TotalTADNum);
    DisListAll = cell(TotalTADNum, TotalTADNum);
    Median = zeros(TotalTADNum, TotalTADNum);
    
    % Initialize counter for progress bar
    steps = TotalTADNum * TotalTADNum;
    varcount = 0;
    for i = 1:TotalTADNum
        for j = 1:TotalTADNum
            varcount = varcount + 1;
            % Update progress bar
            waitbar(varcount/steps, wait_bar, ...
                sprintf('percentage: %2.4f %%', 100 * varcount/steps))
            if getappdata(wait_bar, 'canceling')
                sprintf('User cancelled...exiting...')
                delete(wait_bar)
                return
            end
            % Find all instances of Chr where Chr.r(i) = Chr.r(j) = 1
            tmp = cell2mat(arrayfun(@(x) x.r(i) == 1 & x.r(j) == 1, Chr, 'UniformOutput', false));
            % Compute the distance of Chr instances where the tmp variable is 1
            DisList = arrayfun(@(x) distance(x, i, j), Chr(tmp == 1));
    
            Mean(i, j) = mean(DisList);
            Std(i, j) = std(DisList);
            SEM(i, j) = std(DisList) / (length(DisList))^0.5;
            NofData(i, j) = length(DisList);
            DisListAll{i, j} = DisList;
            Median(i, j) = median(DisList);
        end
    end
    % delete the progress bar and stop the timer
    delete(wait_bar)
    toc
    
    %% make figures
    tic
    
    figure2 = figure;
    colormap(hot);
    % Create axes
    axes1 = axes('Parent', figure2);
    hold(axes1, 'on');
    % Create image
    image(Mean, 'Parent', axes1, 'CDataMapping', 'scaled');
    title(['Mean pairwise distance (um) -'  AgeRange]);
    % Uncomment the following line to preserve the X-limits of the axes
    xlim(axes1, [0.5 40.5]); %28.5
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim(axes1, [0.5 40.5]);
    box(axes1, 'on');
    axis(axes1, 'ij');
    yticks(1:TotalTADNum)
    yticklabels({'16.7', '16.8', '16.9', '17.0', '17.2', '17.4', '17.5', '17.6', '17.7', '17.8', '17.9', '18.0', '18.1', '18.2', '18.3', '18.4', '18.5', '18.6', '18.7', '18.8', '18.9', '19.0', '19.1', '19.2', '19.3', '19.4', '19.5', '19.6', '19.7', '19.8', '19.9', '20.0', '20.1', '20.2', '20.3', '20.4', '20.5', '20.6', '20.7', '20.8', '20.9', '21.0', '21.05'})
    xticks(1:TotalTADNum)
    xticklabels({'16.7', '16.8', '16.9', '17.0', '17.2', '17.4', '17.5', '17.6', '17.7', '17.8', '17.9', '18.0', '18.1', '18.2', '18.3', '18.4', '18.5', '18.6', '18.7', '18.8', '18.9', '19.0', '19.1', '19.2', '19.3', '19.4', '19.5', '19.6', '19.7', '19.8', '19.9', '20.0', '20.1', '20.2', '20.3', '20.4', '20.5', '20.6', '20.7', '20.8', '20.9', '21.0', '21.05'})
    xtickangle(90)
    colorbar('peer', axes1);
    caxis([0 1.3])
    %PlotProp
    axis square
    print('Mean_early','-dsvg')
    if save_flag == 1
        saveas(gcf,[dirName '/Mean pairwise distance (um) -'  AgeRange])
    end
    
    
    figure3 = figure;
    colormap(hot);
    % Create axes
    axes1 = axes('Parent', figure3);
    hold(axes1, 'on');
    % Create image
    image(Median, 'Parent', axes1, 'CDataMapping', 'scaled');
    title(['Median pairwise distance (um) -'  AgeRange]);
    % Uncomment the following line to preserve the X-limits of the axes
    xlim(axes1, [0.5 43.5]);
    % Uncomment the following line to preserve the Y-limits of the axes
    ylim(axes1, [0.5 43.5]);
    box(axes1, 'on');
    axis(axes1, 'ij');
    yticks(1:TotalTADNum)
    yticklabels({'16.7', '16.8', '16.9', '17.0', '17.2', '17.4', '17.5', '17.6', '17.7', '17.8', '17.9', '18.0', '18.1', '18.2', '18.3', '18.4', '18.5', '18.6', '18.7', '18.8', '18.9', '19.0', '19.1', '19.2', '19.3', '19.4', '19.5', '19.6', '19.7', '19.8', '19.9', '20.0', '20.1', '20.2', '20.3', '20.4', '20.5', '20.6', '20.7', '20.8', '20.9', '21.0', '21.05'})
    xticks(1:TotalTADNum)
    xticklabels({'16.7', '16.8', '16.9', '17.0', '17.2', '17.4', '17.5', '17.6', '17.7', '17.8', '17.9', '18.0', '18.1', '18.2', '18.3', '18.4', '18.5', '18.6', '18.7', '18.8', '18.9', '19.0', '19.1', '19.2', '19.3', '19.4', '19.5', '19.6', '19.7', '19.8', '19.9', '20.0', '20.1', '20.2', '20.3', '20.4', '20.5', '20.6', '20.7', '20.8', '20.9', '21.0', '21.05'})
    xtickangle(90)
    colorbar('peer', axes1);
    caxis([0 1.3])
    %PlotProp
    axis square
    print('Median_early','-dsvg')
    if save_flag == 1
        saveas(gcf, [dirName '/Median pairwise distance (um) -'  AgeRange])
    end
    
    figure
    imagesc(SEM)
    colorbar
    title(['SEM of pairwise distance - ' AgeRange]);
    %PlotProp
    axis square
    yticks(1:TotalTADNum)
    yticklabels({'16.7', '16.8', '16.9', '17.0', '17.2', '17.4', '17.5', '17.6', '17.7', '17.8', '17.9', '18.0', '18.1', '18.2', '18.3', '18.4', '18.5', '18.6', '18.7', '18.8', '18.9', '19.0', '19.1', '19.2', '19.3', '19.4', '19.5', '19.6', '19.7', '19.8', '19.9', '20.0', '20.1', '20.2', '20.3', '20.4', '20.5', '20.6', '20.7', '20.8', '20.9', '21.0', '21.05'})
    xticks(1:TotalTADNum)
    xticklabels({'16.7', '16.8', '16.9', '17.0', '17.2', '17.4', '17.5', '17.6', '17.7', '17.8', '17.9', '18.0', '18.1', '18.2', '18.3', '18.4', '18.5', '18.6', '18.7', '18.8', '18.9', '19.0', '19.1', '19.2', '19.3', '19.4', '19.5', '19.6', '19.7', '19.8', '19.9', '20.0', '20.1', '20.2', '20.3', '20.4', '20.5', '20.6', '20.7', '20.8', '20.9', '21.0', '21.05'})
    xtickangle(90)
    print('SEM_early','-dsvg')
    if save_flag == 1
        saveas(gcf, [dirName '/SEM of pairwise distance - ' AgeRange])
    end
    
    
    figure
    imagesc(Std)
    colorbar
    title(['Std of pairwise distance - ' AgeRange]);
    %PlotProp
    axis square
    yticks(1:TotalTADNum)
    yticklabels({'16.7', '16.8', '16.9', '17.0', '17.2', '17.4', '17.5', '17.6', '17.7', '17.8', '17.9', '18.0', '18.1', '18.2', '18.3', '18.4', '18.5', '18.6', '18.7', '18.8', '18.9', '19.0', '19.1', '19.2', '19.3', '19.4', '19.5', '19.6', '19.7', '19.8', '19.9', '20.0', '20.1', '20.2', '20.3', '20.4', '20.5', '20.6', '20.7', '20.8', '20.9', '21.0', '21.05'})
    xticks(1:TotalTADNum)
    xticklabels({'16.7', '16.8', '16.9', '17.0', '17.2', '17.4', '17.5', '17.6', '17.7', '17.8', '17.9', '18.0', '18.1', '18.2', '18.3', '18.4', '18.5', '18.6', '18.7', '18.8', '18.9', '19.0', '19.1', '19.2', '19.3', '19.4', '19.5', '19.6', '19.7', '19.8', '19.9', '20.0', '20.1', '20.2', '20.3', '20.4', '20.5', '20.6', '20.7', '20.8', '20.9', '21.0', '21.05'})
    xtickangle(90)
    print('STD_early','-dsvg')
    if save_flag == 1
        saveas(gcf, [dirName '/Std of pairwise distance - ' AgeRange])
    end
    
    
    figure
    imagesc(NofData)
    colorbar
    title(['Number of measurements - ' AgeRange]);
    %PlotProp
    axis square
    print('NofDATA_early','-dsvg')
    if save_flag == 1
        saveas(gcf, [dirName '/Number of measurements - ' AgeRange])
    end
    
    
    histogram([Chr.TracesInTerritory])
    title('wild-type segmentation results')
    ylabel('number of traces')
    xlabel('total traces in territory volume')
    %PlotProp
    if save_flag == 1
        saveas(gcf, [dirName '/wild-type traces in terry-segmentation results' AgeRange])
    end
    
    
    TracesInTerritory = [Chr.TracesInTerritory].';
    TracesPerSegment = [];
    TracesPerSegment(1, 1) = sum(TracesInTerritory(:) == 1) / length(Chr) * 100;
    TracesPerSegment(1, 2) = sum(TracesInTerritory(:) == 2) / length(Chr) * 100;
    TracesPerSegment(1, 3) = sum(TracesInTerritory(:) == 3) / length(Chr) * 100;
    TracesPerSegment(1, 4) = sum(TracesInTerritory(:) == 4) / length(Chr) * 100;
    figure
    bar(TracesPerSegment)
    title('wild-type segmentation results')
    ylabel('% of total traces')
    xlabel('traces within territory volume')
    ylim([0 100])
    %PlotProp
    print('wt_segmentation_early','-dsvg')
    if save_flag == 1
        saveas(gcf, [dirName '/wild-type segmentation results' AgeRange])
    end
end
toc

%% Save matrices
tic

if save_flag == 1
    disp('Saving Mean Matrix...')
    save([dirName '/MeanPairDistanceMatrix' ...
        num2str(SegChannel) '.mat'], 'Mean');
    disp('Saving Median Matrix...')
    save([dirName '/MedianPairDistanceMatrix' ...
        num2str(SegChannel) '.mat'], 'Median');
    disp('Saving Standard Deviation Matrix...')
    save([dirName '/StdPairDistanceMatrix' ...
        num2str(SegChannel) '.mat'], 'Std');
    disp('Saving SEM Matrix...')
    save([dirName '/SEMPairDistanceMatrix' ...
        num2str(SegChannel) '.mat'], 'SEM');
    disp('Saving NofData Matrix...')
    save([dirName '/NofDataPairDistanceMatrix' ...
        num2str(SegChannel) '.mat'], 'NofData');
    disp('Saving Chr Matrix...')
    save([dirName '/AllChromosomes' ...
        num2str(SegChannel) '.mat'], 'Chr');
    disp('Saving DisListAll Matrix...')
    save([dirName '/DisListAll' ...
        num2str(SegChannel) '.mat'], 'DisListAll');
    
    % Save L and NucL if we're doing so
    if l_flag == 1
        if SegChannel == 647
            disp('Saving L2 Matrix...')
            save([dirName '/AllL2_' ...
                num2str(SegChannel) '.mat'], 'L2s', '-v7.3');
        else
            disp('Saving L1 Matrix...')
            save([dirName '/AllL1_' ...
                num2str(SegChannel) '.mat'], 'L1s', '-v7.3');
        end
    end
    if n_flag == 1
        disp('Saving NucL Matrix...')
        save([dirName '/AllNucL_' ...
            num2str(SegChannel) '.mat'], 'NucLs', '-v7.3');
    end
end

toc
%return
%% Define functions
function t_or_f = check_prereqs(Trace, TAD_id, index)
    t_or_f = iscell(Trace) && iscell(TAD_id) && length(TAD_id{index}) > ...
        1 && length(TAD_id{1, index}) == length(unique(TAD_id{1, index}));
end

function yy = rog(Trace)
    Mean = mean(Trace, 1); % finds the average coordinate in x, y, z
    Sum = 0;
    N = 0;
    for i = 1:size(Trace, 1) % 1 to number of coordinates in Trace
        Sum = Sum + sum((Trace(i, :) - Mean).^2);
        N = N + 1;
    end
    yy = (Sum / N)^0.5;
end

function dist = distance(vect, i, j)
    xsq = (vect.x(i) - vect.x(j))^2;
    ysq = (vect.y(i) - vect.y(j))^2;
    zsq = (vect.z(i) - vect.z(j))^2;
    dist = sqrt(xsq + ysq + zsq);
end

function flag = check_save_input()
    answer = questdlg(['Would you like to save the matrices to the disk'...
        ' or only to workspace?'], 'Save Matrices', 'Save to Workspace',...
        'Save to Disk', 'Save to Disk');
    switch answer
        case 'Save to Disk'
            disp('Will save to disk')
            flag = 1;
        case 'Save to Workspace'
            disp('Will only save to workspace')
            flag = 0;
    end
end

function check_folder_input()
    answer = questdlg(['The folder already exists, are you sure you'...
        ' would like to continue and overwrite the information?'], ...
        'Warning: Continue overwriting', 'Yes', 'No', 'No');
    switch answer
        case 'Yes'
            disp('Continuing...')
        case 'No'
            error('Exiting, change folder name or check file contents')
    end
end

function flag = check_C_input()
    answer = questdlg('Would you like to save the Chr variable?', ...
        'Chr', 'Yes', 'No', 'No');
    switch answer
        case 'Yes'
            disp('Reading Chr')
            flag = 1;
        case 'No'
            disp('Omiting Chr')
            flag = 0;
    end
end

function flag = check_L_input()
    answer = questdlg('Would you like to save the L2 variable?', ...
        'L2', 'Yes', 'No', 'No');
    switch answer
        case 'Yes'
            disp('Reading L2')
            flag = 1;
        case 'No'
            disp('Omit L2')
            flag = 0;
    end
end

function flag = check_N_input()
    answer = questdlg('Would you like to save the NucL variable?', ...
        'NucL', 'Yes', 'No', 'No');
    switch answer
        case 'Yes'
            disp('Reading NucL')
            flag = 1;
        case 'No'
            disp('Omit NucL')
            flag = 0;
    end
end