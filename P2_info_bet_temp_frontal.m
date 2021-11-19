% calculate information flow between temporal(8 12 13) and frontal (1 3 4)

num_patients = 10;
VOT_conditions =3;
nch=30;
nwin = 43;

temporal_to_frontal = cell(1,num_patients);
frontal_to_temporal = cell(1,num_patients);

temporal_elec = [8 12 13];
frontal_elec = [1 3 4];

for ipat = 1:num_patients % for each patient, run the loop
    VOT_ch_time_tf = zeros(VOT_conditions, nwin); %temporal to frontal
    VOT_ch_time_ft = zeros(VOT_conditions, nwin); %frontal to temporal

    % VOT condition 1 : includes 1,2,3
    cond =1; 
    filename = sprintf('di_matrix_subj%d_VOT123_mem2K4down15ver2.mat', ipat);
    load(filename);
   
    for nt = 1:nwin % calculate degrees for each time window
        mat = di_matrix{nt,1};
        mat(mat<0)=0;
        mat_tf = mat(frontal_elec,temporal_elec); % matrix of info from temp to fron
        mat_ft = mat(temporal_elec,frontal_elec); % matrix of info from frontal to temp

        VOT_ch_time_tf(cond,nt) = sum(sum(mat_tf)); %along the row, in-degrees
        VOT_ch_time_ft(cond,nt) = sum(sum(mat_ft));
    end
    
    % VOT condition 2 : includes 4,5,6
    cond =2; 
    filename = sprintf('di_matrix_subj%d_VOT456_mem2K4down15ver2.mat', ipat);
    load(filename);
   
    for nt = 1:nwin % calculate degrees for each time window
        mat = di_matrix{nt,1};
        mat(mat<0)=0;
        mat_tf = mat(frontal_elec,temporal_elec); % matrix of info from temp to fron
        mat_ft = mat(temporal_elec,frontal_elec); % matrix of info from frontal to temp

        VOT_ch_time_tf(cond,nt) = sum(sum(mat_tf)); %along the row, in-degrees
        VOT_ch_time_ft(cond,nt) = sum(sum(mat_ft));
    end
       
    % VOT condition 3 : includes 7,8,9
    cond =3; 
    filename = sprintf('di_matrix_subj%d_VOT789_mem2K4down15ver2.mat', ipat);
    load(filename);
   
    for nt = 1:nwin % calculate degrees for each time window
        mat = di_matrix{nt,1};
        mat(mat<0)=0;
        mat_tf = mat(frontal_elec,temporal_elec); % matrix of info from temp to fron
        mat_ft = mat(temporal_elec,frontal_elec); % matrix of info from frontal to temp

        VOT_ch_time_tf(cond,nt) = sum(sum(mat_tf)); %along the row, in-degrees
        VOT_ch_time_ft(cond,nt) = sum(sum(mat_ft));
    end
    
    % now save the VOT_ch_time structure to the cell node_indeg
    temporal_to_frontal{1,ipat} = VOT_ch_time_tf;
    frontal_to_temporal{1,ipat} = VOT_ch_time_ft;
end

%% look at temporal to frontal in each time window, averaging across patients
ch_comp_baseline = zeros(num_patients,VOT_conditions, nwin-2);
for ipat = 1:num_patients
    mat = temporal_to_frontal{1,ipat};
    ch_mat = mat(:,3:nwin) - repmat(mat(:,1),1,nwin-2);
    ch_comp_baseline(ipat,:,:) = ch_mat;
end

figure;
for i = 1:nwin-2 % in each time window
    
    in_con1 = squeeze(ch_comp_baseline(:,1,i));%condition 1 (VOT 1 2 3) accross all patients
    in_con2 = squeeze(ch_comp_baseline(:,2,i));%condition 2 (VOT 4 5 6) accross all patients
    in_con3 = squeeze(ch_comp_baseline(:,3,i));%condition 3 (VOT 7 8 9) accross all patients

    errorbar([mean(in_con1) mean(in_con2) mean(in_con3)],[std(in_con1) std(in_con2) std(in_con3)]); 
    waitforbuttonpress
end

%% look at frontal to temporal in each time window
ch_comp_baseline = zeros(num_patients,VOT_conditions, nwin-2);
for ipat = 1:num_patients
    mat = frontal_to_temporal{1,ipat};
    ch_mat = mat(:,3:nwin) - repmat(mat(:,1),1,nwin-2);
    ch_comp_baseline(ipat,:,:) = ch_mat;
end

figure;
for i = 1:nwin-2 % in each time window
    
    in_con1 = squeeze(ch_comp_baseline(:,1,i));%condition 1 (VOT 1 2 3) accross all patients
    in_con2 = squeeze(ch_comp_baseline(:,2,i));%condition 2 (VOT 4 5 6) accross all patients
    in_con3 = squeeze(ch_comp_baseline(:,3,i));%condition 3 (VOT 7 8 9) accross all patients

    errorbar([mean(in_con1) mean(in_con2) mean(in_con3)],[std(in_con1) std(in_con2) std(in_con3)]); 
    waitforbuttonpress
end

%% look at frontal to temporal in each condition, average of each window
ch_comp_baseline = zeros(num_patients,VOT_conditions, nwin-2);
for ipat = 1:num_patients
    mat = frontal_to_temporal{1,ipat};
    ch_mat = mat(:,3:nwin);% - repmat(mat(:,1),1,nwin-2);
    ch_comp_baseline(ipat,:,:) = ch_mat;
end

figure;
for i = 1:num_patients % in each time window
    
    in_con1 = squeeze(ch_comp_baseline(i,1,:));%condition 1 (VOT 1 2 3) accross all patients
    in_con2 = squeeze(ch_comp_baseline(i,2,:));%condition 2 (VOT 4 5 6) accross all patients
    in_con3 = squeeze(ch_comp_baseline(i,3,:));%condition 3 (VOT 7 8 9) accross all patients

    errorbar([mean(in_con1) mean(in_con2) mean(in_con3)],[std(in_con1) std(in_con2) std(in_con3)]); 
    waitforbuttonpress
end
        
%% lets look at the data
figure; 

for ipat = 1:num_patients
    mat = temporal_to_frontal{1,ipat};
    subplot(1,3,1); bar(squeeze(mat(1,:,:))); colorbar; colormap jet;
    subplot(1,3,2); bar(squeeze(mat(2,:,:))); colorbar; colormap jet;
    subplot(1,3,3); bar(squeeze(mat(3,:,:))); colorbar; colormap jet;
    waitforbuttonpress;
end
    
figure; 

for ipat = 1:num_patients
    mat = frontal_to_temporal{1,ipat};
    subplot(1,3,1); bar(squeeze(mat(1,:,:))); colorbar; colormap jet;
    subplot(1,3,2); bar(squeeze(mat(2,:,:))); colorbar; colormap jet;
    subplot(1,3,3); bar(squeeze(mat(3,:,:))); colorbar; colormap jet;
    waitforbuttonpress;
end
    