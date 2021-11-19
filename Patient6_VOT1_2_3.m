%% load the undergrad data
load('21SubjCell_ElecMatrixSeparatedfromTrialInfo.mat');
% this contains subjCell

patient_number = 6; bin1 = 27; bin2 = 28; bin3 = 29;
data1 = extract_VOT_data(subjCell, patient_number,bin1);
data2 = extract_VOT_data(subjCell, patient_number,bin2);
data3 = extract_VOT_data(subjCell, patient_number,bin3);

data = [data1 data2 data3];

%% get data into ch x trials x time format
ntrials = size(data,2);
nch = 30;
ntime = 350;
ch_trials_time_all = zeros(nch,ntrials, ntime);
for i = 1:ntrials
    ch_trials_time_all(:,i,:) = data{1,i};
end

%% lets start with 50ms window sizes
winsize = 50;
nperm =20;

nwins = [1 51 101:5:301]; % these are the start time of each window

count = length(nwins);
di_matrix = cell(count,1);

for index = 1:count % for each window, calculate pairwise DI
tic
time_win =nwins(index) : nwins(index)+ winsize -1;
%extracting data in the relevant time window

ch_trials_time = ch_trials_time_all(:,:,time_win); %extract ch x trials x time in relevant window

connectivity_matrix1 = zeros(nch,nch);
connectivity_matrix2 = zeros(nch,nch);

for fch = 2 : nch %start with from channels - fch
    parfor tch = 1 : fch-1 % to channels start from 1, and end 1 less than fch
        %we want to fill in a matrix, fch is row, tch is column. 
        K=4; %knn parameter
        
        % original DI from x to y and y to x
        x1 = reshape(ch_trials_time(fch,:,:),ntrials,winsize); %from ch
        x2 = reshape(ch_trials_time(tch,:,:),ntrials,winsize); %to channel

       [dixy, diyx]=DIver2(x1,x2,2,K,15);
       
       % DI for bootstrapped data
       temp_xy=0;
       temp_yx=0;
       for nn = 1:nperm
           
           xt2 = zeros(size(x2)); % shuffle time in one time series, for each trial
           %shuffle the time values of signal x2
           for tt = 1:ntrials
               [yy,idx] = datasample(1:winsize,winsize,'Replace', false);
               xt2(tt,:) = x2(tt,idx);
           end
           [dixy_t1, diyx_t1]=DIver2(x1,xt2,2,K,15);
           
           temp_xy = temp_xy + dixy_t1;
           temp_yx = temp_yx + diyx_t1;
       end %end permutations

       dixy = dixy - temp_xy/nperm; % bias correction of DI value using the average permuted DI value
       diyx = diyx - temp_yx/nperm;  
       
       connectivity_matrix1(fch,tch) = dixy;
       connectivity_matrix2(fch,tch) = diyx;
    end%tch
   
end%fch
di_matrix{index,1} = connectivity_matrix1 + connectivity_matrix2';

toc
end% end time windows, index

save('di_matrix_subj6_VOT123_mem2K4down15ver2.mat','di_matrix', '-v7.3'); %with mem 2
