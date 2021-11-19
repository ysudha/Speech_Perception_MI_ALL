% calculate node degrees for all patients

num_patients = 21;
VOT_conditions =3;
nch=30;
nwin = 43;

dim = 2; %2 is indegree, 1 is outdegree
node_indeg = cell(1,num_patients);

for ipat = 1:num_patients % for each patient, run the loop
    VOT_ch_time = zeros(VOT_conditions, nch, nwin);
    
    % VOT condition 1 : includes 1,2,3
    cond =1; 
    filename = sprintf('di_matrix_subj%d_VOT123_mem2K4down15ver2.mat', ipat);
    load(filename);
   
    for nt = 1:nwin % calculate degrees for each time window
        mat = di_matrix{nt,1};
        mat(mat<0)=0;
        VOT_ch_time(cond,:,nt) = sum(mat,dim); %along the row, in-degrees
    end
    
    % VOT condition 2 : includes 4,5,6
    cond =2; 
    filename = sprintf('di_matrix_subj%d_VOT456_mem2K4down15ver2.mat', ipat);
    load(filename);
   
    for nt = 1:nwin % calculate degrees for each time window
        mat = di_matrix{nt,1};
        mat(mat<0)=0;
        VOT_ch_time(cond,:,nt) = sum(mat,dim); %along the row, in-degrees
    end
       
    % VOT condition 3 : includes 7,8,9
    cond =3; 
    filename = sprintf('di_matrix_subj%d_VOT789_mem2K4down15ver2.mat', ipat);
    load(filename);
   
    for nt = 1:nwin % calculate degrees for each time window
        mat = di_matrix{nt,1};
        mat(mat<0)=0;           
        VOT_ch_time(cond,:,nt) = sum(mat,dim); %along the row, in-degrees
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
    

    