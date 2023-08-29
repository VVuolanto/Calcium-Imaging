clear
close all

%% Signal Initialization
list = dir('P:\Calcium Data Analysis\Emre\Data files\PFF+3d\control\*.csv'); %Insert the folder where .csv files of regions are

for n = 1 : size(list) %n regions imaged

    
RawSignal = readmatrix([list(n).folder '\' list(n).name]); %

Background = RawSignal(3:end,end-2); %taking the background fluorescence (found in the third to last column, i.e. last ROI from ImageJ)

CAsignal = RawSignal(3:end,1:end-3)-Background; %Substracting background

for i = 1:size(CAsignal,2) %neurons (columns) detected in the region of interest
    %% 10 point MA filter

    wts = [1/20; repmat(1/10,9,1); 1/20];  
    
    filteredCA = conv(CAsignal(1:end,i),wts,'valid'); %filter Ca signal
    
    %% dF/F0

    %the percentage of signal amplitude below which signal should be
    %considered to be baseline (here 20%); above this, spiking activity 
    %occurs
    fraction_of_max_activity = (max(filteredCA) - min(filteredCA))*0.20; 

    %return only lower values of signal = baseline
    baseline = filteredCA;
    baseline(baseline > min(filteredCA) + fraction_of_max_activity) = min(filteredCA) + fraction_of_max_activity; 
    
    %Calculate DeltaF of F0
    F0 = mean(baseline);
    FF0signal = (filteredCA-F0)/F0;

    %% Normalizing any signal overall trend by substracting a linear regression of the signal
    

    x = 1:length(FF0signal);  
    linregr = [ones(length(x),1) x']\FF0signal;
    Lin_detrend = linregr(2)*x+linregr(1); 
    detrend = FF0signal-Lin_detrend'; %overall detrend the data with a linear approximation
    smooth = smoothdata(detrend,'gaussian',10); %smooth data with a Gaussian filter to improve peak detection
    normalized_smooth = smooth-min(smooth); %translate minimum value of signal to 0


    %% Calculate and store neuron signal peaks

    [CAamps, CAlocs, Peak_Width, Prominence] = findpeaks(normalized_smooth,1,'MinPeakProminence',0.05); %find peaks and their parameters from the filtered signal
  
    if length(CAlocs) > 3 %Exclude signal points that are not active (counts only neurons with more than [qualifier] spikes)

        CA_p2p = diff(CAlocs); %calculate peak-to-peak distances
        CA_p2p = CA_p2p * 0.028; %converting peak-to-peak distance from frames to seconds

        Peak_Width = Peak_Width * 0.028; %converting peak width from frames to seconds

        Prom_Cell{1,i} = mean(Prominence); %storing average prominence of spikes from every neuron in the ROI

        Amp_Cell{1,i} = mean(CAamps); %storing average amplitude of spikes from every neuron in the ROI

        p2p_Cell{1,i} = mean(CA_p2p); %storing average peak-to-peak values from every neuron in the ROI

        Count_Cell{1,i} = length(CAlocs)-2; %storing count of spikes from every neuron in the ROI

        Peak_Width_Cell{1,i} = mean(Peak_Width); %storing average of peak widths from every neuron in the ROI


        
    end
end
%% Export Data

%convert cells into vectors for exporting
Average_Amps = cell2mat(Amp_Cell);
Average_Prom = cell2mat(Prom_Cell);
Average_p2p = cell2mat(p2p_Cell);
Average_Count = cell2mat(Count_Cell);
Average_PW = cell2mat(Peak_Width_Cell);

% Preallocate matrix with all data from a single ROI, capped at 150 neurons detected.
% Add an identifier for every column
ROI_matrix = zeros(151,4); 
ROI_matrix(1,1) = 1001; %AverageAmps = 1001
ROI_matrix(1,2) = 2002; %AverageProm = 2002
ROI_matrix(1,3) = 3003; %AveragePeaks = 3003
ROI_matrix(1,4) = 4004; %PeakCount = 4004
ROI_matrix(1,5) = 5005; %AvgPeakWidth = 5005

% Add data from 'Average' vectors into the single ROI matrix
for k = 1 : length(Average_Amps)
    ROI_matrix(k+1,1) = Average_Amps(k);
    ROI_matrix(k+1,2) = Average_Prom(k);
    ROI_matrix(k+1,3) = Average_p2p(k);
    ROI_matrix(k+1,4) = Average_Count(k);
    ROI_matrix(k+1,5) = Average_PW(k);  
end

% Store each ROI into cells
Group_Data_Cell{n,1} = ROI_matrix;
end

Group_Data_Matrix = cell2mat(Group_Data_Cell);

% Export data into an excel sheet
writematrix(Group_Data_Matrix, '[folder path (with a new filename or an existing filename at the end)]', 'Sheet', '[Which Sheet in Excel you want to paste the data in]', 'Range', '[starting cell in Excel sheet]');
