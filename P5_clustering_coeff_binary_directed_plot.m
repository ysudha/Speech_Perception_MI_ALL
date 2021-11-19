aa = zeros(1,3);
bb = zeros(1,3);
cc = zeros(1,3);
dd = zeros(1,3);
ee = zeros(1,3);
ff = zeros(1,3);
gg = zeros(1,3);
hh = zeros(1,3);

for patient_num=1:21

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

for i = 3:43   %first 2 windows are from the baseline
    temp = (a1.di_matrix{i,1} );% - a3.di_matrix{i,1});caxis([0 0.3]);
    temp(temp< thresh) =0;
   temp=double((temp~=0));
    cluscoef1(:,i) = clustering_coef_bd(temp);
    
    temp = (a2.di_matrix{i,1} );% - a3.di_matrix{i,1});caxis([0 0.3]);
    temp(temp< thresh) =0;
    temp=double((temp~=0));
    cluscoef2(:,i) = clustering_coef_bd(temp);
    
    temp = (a3.di_matrix{i,1} );% - a3.di_matrix{i,1});caxis([0 0.3]);
    temp(temp< thresh) =0;
   temp=double((temp~=0));
    cluscoef3(:,i) = clustering_coef_bd(temp);
    
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

% LT 1:10
mat1 = mean(cluscoef1([1 3 4 8 12 13 7],[1:10]),2);
mat2 = mean(cluscoef2([1 3 4 8 12 13 7],[1:10]),2);
mat3 = mean(cluscoef3([1 3 4 8 12 13 7],[1:10]),2);
[mat1 mat2 mat3]
elec=[4:6];
aa(patient_num,:) = [mean(mat1(elec)) mean(mat2(elec)) mean(mat3(elec))]

% LF 1:10
mat1 = mean(cluscoef1([1 3 4 8 12 13 7],[1:10]),2);
mat2 = mean(cluscoef2([1 3 4 8 12 13 7],[1:10]),2);
mat3 = mean(cluscoef3([1 3 4 8 12 13 7],[1:10]),2);
[mat1 mat2 mat3]
elec=[1:3];
ee(patient_num,:) = [mean(mat1(elec)) mean(mat2(elec)) mean(mat3(elec))]

%LF 11:16
mat1 = mean(cluscoef1([1 3 4 8 12 13 7],[11:16]),2);
mat2 = mean(cluscoef2([1 3 4 8 12 13 7],[11:16]),2);
mat3 = mean(cluscoef3([1 3 4 8 12 13 7],[11:16]),2);
[mat1 mat2 mat3]
elec=[1:3];
cc(patient_num,:) = [mean(mat1(elec)) mean(mat2(elec)) mean(mat3(elec))]

%LT 11:16
mat1 = mean(cluscoef1([1 3 4 8 12 13 7],[11:16]),2);
mat2 = mean(cluscoef2([1 3 4 8 12 13 7],[11:16]),2);
mat3 = mean(cluscoef3([1 3 4 8 12 13 7],[11:16]),2);
[mat1 mat2 mat3]
elec=[4:6];
bb(patient_num,:) = [mean(mat1(elec)) mean(mat2(elec)) mean(mat3(elec))]

%LT 22:27
mat1 = mean(cluscoef1([1 3 4 8 12 13 7],[22:27]),2);
mat2 = mean(cluscoef2([1 3 4 8 12 13 7],[22:27]),2);
mat3 = mean(cluscoef3([1 3 4 8 12 13 7],[22:27]),2);
[mat1 mat2 mat3]
elec=[4:6];
gg(patient_num,:) = [mean(mat1(elec)) mean(mat2(elec)) mean(mat3(elec))]

%LF 22:27
mat1 = mean(cluscoef1([1 3 4 8 12 13 7],[22:27]),2);
mat2 = mean(cluscoef2([1 3 4 8 12 13 7],[22:27]),2);
mat3 = mean(cluscoef3([1 3 4 8 12 13 7],[22:27]),2);
[mat1 mat2 mat3]
elec=[1:3];
hh(patient_num,:) = [mean(mat1(elec)) mean(mat2(elec)) mean(mat3(elec))]


%RT RF 11:16
mat1 = mean(cluscoef1([2 6 7 10 11],[11:16]),2);
mat2 = mean(cluscoef2([2 6 7 10 11],[11:16]),2);
mat3 = mean(cluscoef3([2 6 7 10 11],[11:16]),2);
[mat1 mat2 mat3]
 elec=[1:5];
 dd(patient_num,:) = [mean(mat1(elec)) mean(mat2(elec)) mean(mat3(elec))]

 %RT RF 1:10
mat1 = mean(cluscoef1([2 6 7 10 11],[1:10]),2);
mat2 = mean(cluscoef2([2 6 7 10 11],[1:10]),2);
mat3 = mean(cluscoef3([2 6 7 10 11],[1:10]),2);
[mat1 mat2 mat3]
 elec=[1:5];
 ff(patient_num,:) = [mean(mat1(elec)) mean(mat2(elec)) mean(mat3(elec))]


end

%% plotting
%errorbar(mean(aa),std(aa));
% figure;plot(aa');

fig1=figure; 
%bar(mean(aa));ylim([0,0.5];
aa = normalize_multiple_patients(aa);
barwitherr(std(aa),mean(aa),'FaceColor',[0 0 1]);ylim([0,0.5])
title('Left temporal : 0 to 95 ms'); 
labels_bar={'shortVOTs'; 'midVOTs'; 'longVOTs' };
set(gca,'xticklabel',labels_bar,'FontSize', 20)
str = sprintf('Left_temporal_0_to_95_ms');
box off;
custom_mask = -1*(aa(:,1)) + 2*(aa(:,2)) -1*(aa(:,3));
[h,p] = ttest(custom_mask)
print(fig1,str, '-depsc', '-r0');

fig1=figure;
bb = normalize_multiple_patients(bb);
barwitherr(std(bb),mean(bb),'FaceColor',[0 0 1]);ylim([0,0.5])
title('Left temporal : 50 to 125 ms'); 
labels_bar={'shortVOTs'; 'midVOTs'; 'longVOTs' };
set(gca,'xticklabel',labels_bar,'FontSize', 20)
str = sprintf('Left_temporal_50_to_125_ms');
box off;
custom_mask = -1*(bb(:,1)) + 2*(bb(:,2)) -1*(bb(:,3));
[h,p] = ttest(custom_mask)
print(fig1,str, '-depsc', '-r0');

fig1=figure;
gg = normalize_multiple_patients(gg);
barwitherr(std(gg),mean(gg),'FaceColor',[0 0 1]);ylim([0,0.5])
title('Left temporal : 105 to 180 ms'); 
labels_bar={'shortVOTs'; 'midVOTs'; 'longVOTs' };
set(gca,'xticklabel',labels_bar,'FontSize', 20)
str = sprintf('Left_temporal_105_to_180_ms');
box off;
custom_mask = -1*(gg(:,1)) + 2*(gg(:,2)) -1*(gg(:,3));
[h,p] = ttest(custom_mask)
print(fig1,str, '-depsc', '-r0');

fig1=figure;
hh = normalize_multiple_patients(hh);
barwitherr(std(hh),mean(hh),'FaceColor',[0 0 1]);ylim([0,0.5])
title('Left frontal : 105 to 180 ms'); 
labels_bar={'shortVOTs'; 'midVOTs'; 'longVOTs' };
set(gca,'xticklabel',labels_bar,'FontSize', 20)
str = sprintf('Left_frontal_105_to_180_ms');
box off;
custom_mask = -1*(hh(:,1)) + 2*(hh(:,2)) -1*(hh(:,3));
[h,p] = ttest(custom_mask)
print(fig1,str, '-depsc', '-r0');


fig1 = figure;
cc = normalize_multiple_patients(cc);
barwitherr(std(cc),mean(cc),'FaceColor',[1 0 0]);ylim([0,0.5])
title('Left frontal : 50 to 125 ms'); 
labels_bar={'shortVOTs'; 'midVOTs'; 'longVOTs' };
set(gca,'xticklabel',labels_bar,'FontSize', 20)
str = sprintf('Left_frontal_50_to_125_ms');
box off;
custom_mask = -1*(cc(:,1)) + 2*(cc(:,2)) -1*(cc(:,3));
[h,p] = ttest(custom_mask)
print(fig1,str, '-depsc', '-r0');

fig1 = figure;
ee = normalize_multiple_patients(ee);
barwitherr(std(ee),mean(ee),'FaceColor',[1 0 0]);ylim([0,0.5])
title('Left frontal : 0 to 95 ms'); 
labels_bar={'shortVOTs'; 'midVOTs'; 'longVOTs' };
set(gca,'xticklabel',labels_bar,'FontSize', 20)
str = sprintf('Left_frontal_0_to_95_ms');
box off;
custom_mask = -1*(ee(:,1)) + 2*(ee(:,2)) -1*(ee(:,3));
[h,p] = ttest(custom_mask)
print(fig1,str, '-depsc', '-r0');

fig1 = figure;
dd = normalize_multiple_patients(dd);
barwitherr(std(dd),mean(dd),'FaceColor',[1 0 1]);ylim([0,0.5])
title('Right frontal and temporal : 50 to 125 ms'); 
labels_bar={'shortVOTs'; 'midVOTs'; 'longVOTs' };
set(gca,'xticklabel',labels_bar,'FontSize', 20)
str = sprintf('Right_frontal_temporal_50_to_125_ms');
box off;
custom_mask = -1*(dd(:,1)) + 2*(dd(:,2)) -1*(dd(:,3));
[h,p] = ttest(custom_mask)
print(fig1,str, '-depsc', '-r0');

fig1 = figure;
ff = normalize_multiple_patients(ff);
barwitherr(std(ff),mean(ff),'FaceColor',[1 0 1]);ylim([0,0.5])
title('Right frontal and temporal : 0 to 95 ms'); 
labels_bar={'shortVOTs'; 'midVOTs'; 'longVOTs' };
set(gca,'xticklabel',labels_bar,'FontSize', 20)
str = sprintf('Right_frontal_temporal_0_to_95_ms');
box off;
custom_mask = -1*(ff(:,1)) + 2*(ff(:,2)) -1*(ff(:,3));
[h,p] = ttest(custom_mask)
print(fig1,str, '-depsc', '-r0');
