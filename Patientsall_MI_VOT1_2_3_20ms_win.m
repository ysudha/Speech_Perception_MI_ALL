%% load the undergrad data
load('21SubjCell_ElecMatrixSeparatedfromTrialInfo.mat');
% this contains subjCell

mi_matrix_VOT123_K4 = cell(22,1);
for patient_number = 1:22
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

    %% lets start with 25ms window sizes
    winsize = 10;
    nperm =5;

    nwins = [[1:10:91] [101:10:341]]; % these are the start time of each window

    count = length(nwins);
    mi_matrix = cell(count,1);

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

           [mixy]=MutualInformation(x1,x2,K);

           % DI for bootstrapped data
           temp_xy=0;
           for nn = 1:nperm

               xt2 = zeros(size(x2)); % shuffle time in one time series, for each trial
               %shuffle the time values of signal x2
               for tt = 1:ntrials
                   [yy,idx] = datasample(1:winsize,winsize,'Replace', false);
                   xt2(tt,:) = x2(tt,idx);
               end
               [mixy_t1]=MutualInformation(x1,xt2,K);

               temp_xy = temp_xy + mixy_t1;
           end %end permutations

           mixy = mixy - temp_xy/nperm; % bias correction of MI value using the average permuted DI value


           connectivity_matrix1(fch,tch) = mixy;
           connectivity_matrix2(fch,tch) = mixy;
        end%tch

    end%fch
    mi_matrix{index,1} = connectivity_matrix1 + connectivity_matrix2';

    toc
    end% end time windows, index

    mi_matrix_VOT123_K4{patient_number,1}= mi_matrix;
    save('mi_matrix_VOT123_20ms_win.mat','mi_matrix_VOT123_K4', '-v7.3'); %with mem 2

end

save('mi_matrix_VOT123_20ms_win.mat','mi_matrix_VOT123_K4', '-v7.3'); %with mem 2
