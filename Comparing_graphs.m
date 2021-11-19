%% comparing MI and DI graphs for undergrad 1
% comparing MI and DI graphs to understand the differences, for undergrad 1,
% VOT 123

% load MI data 
% contains cell structure with 2 undergrads (ug), eventually will
% have all 21 undergrads
load('mi_matrix_VOT123_K4.mat');
mi_ug1 = mi_matrix_VOT123_K4{1,1}; % this is the 43 time windows of ug 1

% load DI data
load('di_matrix_subj1_VOT123_mem2K4down15ver2');

for i = 1:43
    figure;
    subplot(1,2,1); imagesc(mi_ug1{i,1}); colormap jet; colorbar; axis square;
    subplot(1,2,2); imagesc(di_matrix{i,1}); colormap jet; colorbar; axis square;
end


%% comparing apahsia and healthy undergrad graphs

% load DI data
ug1=load('di_matrix_subj1_VOT123_mem2K4down15ver2');
%load aphasia data
aph = load('di_matrix_Aphasia_VOT123_mem2K4down15ver2');

% for i = 1:43
%     figure;
%     subplot(1,2,1); imagesc(ug1.di_matrix{i,1}); colormap jet; colorbar; axis square;
%     subplot(1,2,2); imagesc(aph.di_matrix{i,1}); colormap jet; colorbar; axis square;
% end

%% comparing degres of graphs

% calculate node degrees for all patients

num_patients = 21;
VOT_conditions =3;
nch=30;
nwin = 43;

dim = 2; %2 is indegree, 1 is outdegree
node_indeg = cell(1,num_patients+1);

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

% aphasia patient
% VOT condition 1 : includes 1,2,3
cond =1; 
filename = sprintf('di_matrix_Aphasia_VOT123_mem2K4down15ver2.mat');
load(filename);

for nt = 1:nwin % calculate degrees for each time window
    mat = di_matrix{nt,1};
    mat(mat<0)=0;
    VOT_ch_time(cond,:,nt) = sum(mat,dim); %along the row, in-degrees
end

node_indeg{1,num_patients+1} = VOT_ch_time; % last one is aphasia patient
        
% lets look at the data

% for ipat = 1:num_patients+1
%     figure; 
%     mat = node_indeg{1,ipat};
%     imagesc(squeeze(mat(1,:,:))); colorbar; colormap jet;
%     %subplot(1,3,1); imagesc(squeeze(mat(1,:,:))); colorbar; colormap jet;
%     %subplot(1,3,2); imagesc(squeeze(mat(2,:,:))); colorbar; colormap jet;
%     %subplot(1,3,3); imagesc(squeeze(mat(3,:,:))); colorbar; colormap jet;
%     %waitforbuttonpress;
% end

% lets look at in_degrees and out_degrees "expectation" curves
for ipat = 1:num_patients+1
    mat = node_indeg{1,ipat};
    %imagesc(squeeze(mat(1,:,:))); colorbar; colormap jet;
     test = squeeze(mat(1,:,:));
     figure; hold on;
     for j = 1:30
         if all(test(j,1)> test(j,3:43)) || all(test(j,2)> test(j,3:43))
             plot(1:43,test(j,:)); legend j
         end
     end
     hold off;
    
end

%% normalize matrices

% extract all DI matrices in a structure : patient number x ch x ch x
% time

all_DI = zeros(num_patients+1,nch,nch,nwin);

for ipat = 1:num_patients % for each patient, run the loop
   
    % VOT condition 1 : includes 1,2,3
    filename = sprintf('di_matrix_subj%d_VOT123_mem2K4down15ver2.mat', ipat);
    load(filename);
   
    for nt = 1:nwin % calculate degrees for each time window
        mat = di_matrix{nt,1};
        mat(mat<0)=0;
        all_DI(ipat,:,:,nt) = mat; %along the row, in-degrees
    end
  
end

% aphasia patient
% VOT condition 1 : includes 1,2,3
cond =1; 
filename = sprintf('di_matrix_Aphasia_VOT123_mem2K4down15ver2.mat');
load(filename);

for nt = 1:nwin % calculate degrees for each time window
    mat = di_matrix{nt,1};
    mat(mat<0)=0;
    all_DI(num_patients+1,:,:,nt) = mat; %along the row, in-degrees
end

% normalize relative to baseline

for ipat = 1:num_patients+1
    baseline = squeeze(all_DI(ipat,:,:,1));
    all_DI(ipat,:,:,:) = squeeze(all_DI(ipat,:,:,:)) - repmat(baseline,1,1,nwin);
end

% threshold the graphs, based on 2nd baseline window
for ipat = 1:num_patients+1
    person_DI_values = squeeze(all_DI(ipat,:,:,:));
    [thresh,~] = find_threshold_connections(person_DI_values, 2, 5);
    person_DI_values(person_DI_values<thresh)=0;
    all_DI(ipat,:,:,:) = person_DI_values;
end

% calculate in-degrees
in_degrees_pat = zeros(num_patients+1,nch,nwin);
for ipat = 1:num_patients+1
    person_DI_values = squeeze(all_DI(ipat,:,:,:));
    for i = 1:nwin
        in_degrees_pat(ipat,:,i) = sum(squeeze(person_DI_values(:,:,i)),2);
    end
end

% plot indegrees
for ipat = 1:num_patients+1
    in_degrees = squeeze(in_degrees_pat(ipat,:,:));
    figure; imagesc(in_degrees); colormap jet; colorbar;
end

% plot adjacency matrices
for ipat = 1:num_patients+1
    person_DI_values = squeeze(all_DI(ipat,:,:,3));
    figure; imagesc(person_DI_values); colormap jet; colorbar; axis square;
end
