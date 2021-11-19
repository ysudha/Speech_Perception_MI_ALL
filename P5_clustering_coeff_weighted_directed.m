aa = zeros(1,3);
for patient_num=1:10

str1 = sprintf('di_matrix_subj%d_VOT123_mem2K4down15ver2.mat',patient_num);
str2 = sprintf('di_matrix_subj%d_VOT456_mem2K4down15ver2.mat',patient_num);
str3 = sprintf('di_matrix_subj%d_VOT789_mem2K4down15ver2.mat',patient_num);

a1 = load(str1);
a2 = load(str2);
a3 = load(str3);
load('ChanNamesandLocs.mat');
chnames = chanNamesLocs(:,1);

% % lets pick a threshold such such baseline has 5% connection density
a = a1.di_matrix{3,1}; top_per = 5;
b = a(:);
b = sort(b,'descend');
n = ceil(top_per*length(b)/100); % top 3 percent of connections threshold
thresh = b(n);

time_win = [101:5:301];
cluscoef1 = zeros(30,41);
cluscoef2 = zeros(30,41);
cluscoef3 = zeros(30,41);

for i = 3:43   
    temp = (a1.di_matrix{i,1} );% - a3.di_matrix{i,1});caxis([0 0.3]);
    temp(temp< thresh) =0;
    temp = temp./max(temp(:));
   % temp=double((temp~=0));
    cluscoef1(:,i) = clustering_coef_wd(temp);
    
    temp = (a2.di_matrix{i,1} );% - a3.di_matrix{i,1});caxis([0 0.3]);
    temp(temp< thresh) =0;
    temp = temp./max(temp(:));
    %temp=double((temp~=0));
    cluscoef2(:,i) = clustering_coef_wd(temp);
    
    temp = (a3.di_matrix{i,1} );% - a3.di_matrix{i,1});caxis([0 0.3]);
    temp(temp< thresh) =0;
    temp = temp./max(temp(:));
   % temp=double((temp~=0));
    cluscoef3(:,i) = clustering_coef_wd(temp);
    
end

% figure;
% subplot(1,3,1); imagesc(cluscoef1); colormap jet; colorbar;
% subplot(1,3,2); imagesc(cluscoef2); colormap jet; colorbar;
% subplot(1,3,3); imagesc(cluscoef3); colormap jet; colorbar;
% 
% figure;
% subplot(1,3,1); imagesc(cluscoef1([1 3 4 8 12 13],:)); colormap jet; colorbar;
% subplot(1,3,2); imagesc(cluscoef2([1 3 4 8 12 13],:)); colormap jet; colorbar;
% subplot(1,3,3); imagesc(cluscoef3([1 3 4 8 12 13],:)); colormap jet; colorbar;

% lets compare clustering coefficient of frontal electrodes from time
% windows 1 to 10, 11:16, 11 to rest, for each of the conditions

mat1 = sum(cluscoef1([1 3 4 8 12 13],[1:10]),2);
mat2 = sum(cluscoef2([1 3 4 8 12 13],[1:10]),2);
mat3 = sum(cluscoef3([1 3 4 8 12 13],[1:10]),2);
[mat1 mat2 mat3]

mat1 = sum(cluscoef1([1 3 4 8 12 13],[11:16]),2);
mat2 = sum(cluscoef2([1 3 4 8 12 13],[11:16]),2);
mat3 = sum(cluscoef3([1 3 4 8 12 13],[11:16]),2);
[mat1 mat2 mat3]

mat1 = sum(cluscoef1([1 3 4 8 12 13],[11:25]),2);
mat2 = sum(cluscoef2([1 3 4 8 12 13],[11:25]),2);
mat3 = sum(cluscoef3([1 3 4 8 12 13],[11:25]),2);
[mat1 mat2 mat3]
aa(patient_num,:) = [mat1(4,:) mat2(4,:) mat3(4,:)]

end

errorbar(mean(aa),std(aa))
