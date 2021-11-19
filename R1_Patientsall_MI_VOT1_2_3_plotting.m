%% load the data
load('mi_matrix_VOT123_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
for ipat = 1:np
    mi_mat = mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    baseline = mi_mat{1,1}; %first window is the baseline
    mi_ch_ch_time = zeros(nch,nch,nwin);
    mi_concat = []; % concatenating the 43 time windows for plotting
    for iwin = 1:nwin
        mi_ch_ch_time(:,:,iwin) = mi_mat{iwin,1} - baseline;
        mi_concat = [mi_concat (mi_mat{iwin,1} - baseline)];
    end
    figure('Position',[100 1000 2000 500]);
    %mi_concat(mi_concat<0)=0;
    imagesc(mi_concat); colormap jet; colorbar;
end

%% plot raw MI

load('mi_matrix_VOT123_K4.mat'); %with mem 2
np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
for ipat = 1:np
    mi_mat = mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    for iwin = 1:nwin
        figure; imagesc(squeeze(mi_mat{iwin,1})); colormap jet; axis square;
        set(gcf,'color','w');
    end
end
        