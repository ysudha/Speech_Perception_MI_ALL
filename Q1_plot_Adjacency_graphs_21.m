% Program written on Jan 29, 2020, to look at the adjacency matrices of all
% patients in one time window, 2 different thresholdings.

figure; 
for patient_num = 1:21

str1 = sprintf('di_matrix_subj%d_VOT123_mem2K4down15ver2.mat',patient_num);
str2 = sprintf('di_matrix_subj%d_VOT456_mem2K4down15ver2.mat',patient_num);
str3 = sprintf('di_matrix_subj%d_VOT789_mem2K4down15ver2.mat',patient_num);

a1 = load(str1);
a2 = load(str2);
a3 = load(str3);
load('ChanNamesandLocs.mat');
chnames = chanNamesLocs(:,1);

imagesc(squeeze(a1.di_matrix{1,1})); colormap jet; colorbar;
    waitforbuttonpress;
end
