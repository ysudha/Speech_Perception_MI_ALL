%% load the data
miVOT123 = load('mi_matrix_VOT123_K4.mat'); %with mem 2
miVOT456 = load('mi_matrix_VOT456_K4.mat'); 

np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
% plotting the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases
for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    figure('Position',[100 1000 2000 500]);
    for iwin = 1:nwin
       subplot(4,11,iwin); 
       mat = squeeze(mi456{iwin,1}) - squeeze(mi123{iwin,1});
       mat(mat<0)=0;
     % mat = mat([1 3 4 8 12 13],[1 3 4 8 12 13]);
       imagesc(mat); colormap jet; axis square; %colorbar; %caxis([-0.15 0.15]);
       %str = sprintf('%d to %d ms', ntim(iwin)-1, ntim(iwin)+50-1);
       %title(str);
       if iwin ==nwin
           %subplot(4,11,iwin+1); box off; axis off; str = sprintf('Person %d',ipat); title(str); colormap jet;% caxis([-0.15 0.15]);colorbar;       
       end
    end
end

% plotting number of edges that show an increase
ncolors= 22;
redc = [1, 0, 0];
pinkc = [255, 192, 203]/255;
colors_p = [linspace(redc(1),pinkc(1),ncolors)', linspace(redc(2),pinkc(2),ncolors)', linspace(redc(3),pinkc(3),ncolors)'];

time_number_of_edges = zeros(np,nwin);
average_increase = zeros(np,nwin);
    figure('Position',[100 1000 2000 500]);
for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};


    for iwin = 1:nwin
       mat = squeeze(mi456{iwin,1}) - squeeze(mi123{iwin,1});
       mat(mat<0)=0;
       mat = mat([1 3 4 8 12 13],[1 3 4 8 12 13]);
       time_number_of_edges(ipat,iwin) = nnz(mat);
       average_increase(ipat,iwin) = sum(sum(mat))/nnz(mat);
    end
    subplot(1,2,1); plot(1:nwin,time_number_of_edges(ipat,:), 'LineWidth',6, 'color', colors_p(ipat,:)); hold on;
    subplot(1,2,2); plot(1:nwin,average_increase(ipat,:), 'LineWidth',6,'color',colors_p(ipat,:)); hold on;

    str = sprintf('Person %d',ipat); title(str);  
       
end

%time_number_of_edges has info about # edges in the interesting regions

%% creating new graphs, edge = midVOT - lowVOT ...looking at increases only

np = 22; %number of people 21 undergrads, 1 aphasia
nwin = 43;
nch = 30;
ntim = [1 51 101:5:301];
% creating the adjacency matrices, with edge = VOT456- VOT123, only
% positive increases diff_mid_lowVOT
diff_mid_lowVOT = zeros(np,nch,nch,nwin);
for ipat = 1:np
    mi123 = miVOT123.mi_matrix_VOT123_K4{ipat,1}; %contains 43 time windows
    mi456 = miVOT456.mi_matrix_VOT456_K4{ipat,1};
    for iwin = 1:nwin
       mat = squeeze(mi456{iwin,1}) - squeeze(mi123{iwin,1});
       mat(mat<0)=0;
       diff_mid_lowVOT(ipat,:,:,iwin)= mat;
       
    end
end

save('diff_mid_lowVOT','diff_mid_lowVOT','-v7.3');