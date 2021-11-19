% calculate node degrees for all patients

num_patients = 10;
VOT_conditions =3;
nch=30;
nwin = 43;

node_indeg = cell(1,num_patients);

for ipat = 1:num_patients % for each patient, run the loop
    VOT_ch_time = zeros(VOT_conditions, nch, nwin);
    
    % baseline of patient i, window1, VOT123 , making max value 1
    thresh = load('di_matrix_subj1_VOT123_mem2K4down15ver2.mat');
    [normalizing_val ~] = max(max(thresh.di_matrix{1,1}));
    
    % VOT condition 1 : includes 1,2,3
    cond =1; 
    filename = sprintf('di_matrix_subj%d_VOT123_mem2K4down15ver2.mat', ipat);
    load(filename);
   
    for nt = 1:nwin % calculate degrees for each time window
        mat = di_matrix{nt,1}./normalizing_val;
        mat(mat<0)=0;
        VOT_ch_time(cond,:,nt) = sum(mat,2); %along the row, in-degrees
    end
    
    % VOT condition 2 : includes 4,5,6
    cond =2; 
    filename = sprintf('di_matrix_subj%d_VOT456_mem2K4down15ver2.mat', ipat);
    load(filename);
   
    for nt = 1:nwin % calculate degrees for each time window
        mat = di_matrix{nt,1}./normalizing_val;
        mat(mat<0)=0;
        VOT_ch_time(cond,:,nt) = sum(mat,2); %along the row, in-degrees
    end
       
    % VOT condition 3 : includes 7,8,9
    cond =3; 
    filename = sprintf('di_matrix_subj%d_VOT789_mem2K4down15ver2.mat', ipat);
    load(filename);
   
    for nt = 1:nwin % calculate degrees for each time window
        mat = di_matrix{nt,1}./normalizing_val;
        mat(mat<0)=0;           
        VOT_ch_time(cond,:,nt) = sum(mat,2); %along the row, in-degrees
    end
    
    % now save the VOT_ch_time structure to the cell node_indeg
    node_indeg{1,ipat} = VOT_ch_time;
end
        
% lets look at the data
figure; 

for ipat = 1:num_patients
    mat = node_indeg{1,ipat};
    subplot(1,3,1); imagesc(squeeze(mat(1,:,:))); colorbar; colormap jet;
    subplot(1,3,2); imagesc(squeeze(mat(2,:,:))); colorbar; colormap jet;
    subplot(1,3,3); imagesc(squeeze(mat(3,:,:))); colorbar; colormap jet;
    waitforbuttonpress;
end

%% temporal electrodes
figure; 
for ipat = 1:num_patients
    mat = node_indeg{1,ipat};
    elec12 = squeeze(mat(1,12,:));
    elec13 = squeeze(mat(1,13,:));
    elec8 = squeeze(mat(1,8,:));
    hold on;
    plot(elec12,'b-');
    plot(elec13,'r');
    plot(elec8,'g-');
    hold off;
    waitforbuttonpress;
    clf;
end

%% frontal electrodes

time_win = [101-100:5:301-100];

figure; 
for ipat = 1:num_patients
    mat = node_indeg{1,ipat};
    elec1 = squeeze(mat(1,1,:));
    elec3 = squeeze(mat(1,3,:));
    elec4 = squeeze(mat(1,4,:));
    hold on;
    plot([time_win], elec1,'b-');
    plot([time_win],elec3,'r');
    plot([time_win],elec4,'g-');
    hold off;
    waitforbuttonpress;
    clf;
end

%%
figure; ipat=1;
for ch = 1:30
    mat = node_indeg{1,ipat};
    elec = squeeze(mat(1,ch,:));

    plot([time_win],elec(3:end)-elec(1),'r*-');
    waitforbuttonpress;
end
