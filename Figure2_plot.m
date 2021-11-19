% plotting some example EEG signals to make the plot with sliding windows
% used for Figure 2 as well
%% load the undergrad data
load('21SubjCell_ElecMatrixSeparatedfromTrialInfo.mat');
% this contains subjCell

patient_number = 1%:22
if patient_number<22
    bin1 = 27; bin2 = 28; bin3 = 29;
    data1 = extract_VOT_data(subjCell, patient_number,bin1);
    data2 = extract_VOT_data(subjCell, patient_number,bin2);
    data3 = extract_VOT_data(subjCell, patient_number,bin3);
    
    data = [data1 data2 data3];
else
    load('FWEEGTrialsv2.mat');
    bin1 = 5; bin2 = 10; bin3 = 15;
    data1 = extract_VOT_data_aphasia(FWMatrNoNeutral,bin1);
    data2 = extract_VOT_data_aphasia(FWMatrNoNeutral,bin2);
    data3 = extract_VOT_data_aphasia(FWMatrNoNeutral,bin3);
    
    data = [data1 data2 data3];
end

%% get data into ch x trials x time format
ntrials = size(data,2);
nch = 30;
ntime = 350;
ch_trials_time_all = zeros(nch,ntrials, ntime);
for i = 1:ntrials
    ch_trials_time_all(:,i,:) = data{1,i};
end


ch = 6; %choose channel number : 6

trials_time = squeeze(ch_trials_time_all(ch,:,:));

nwins = [1 51 101:5:301]; % these are the start time of each window

%% plotting the figure
figure; hold on;
for i = 2:5
    plot(trials_time(i,:)'+1.5*i*30,'LineWidth',6);
end
    plot(trials_time(1,:)'+0*i*30,'LineWidth',6);

box off;
set(gcf,'color','w');

set(gca,'ytick',[])


