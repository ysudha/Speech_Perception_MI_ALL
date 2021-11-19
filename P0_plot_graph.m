patient_num = 1;

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
thresh = b(n)

figure; imagesc(squeeze(a3.di_matrix{5,1})); colormap jet; colorbar;
%%thresh = 0.1;
figure;
time_win = [101:5:301];
for i = 3:43   
    temp = (a3.di_matrix{i,1});% - a3.di_matrix{i,1});caxis([0 0.3]);
    temp(temp< thresh) =0;
    G = digraph(temp,chnames);
    plot(G,'Layout','circle','LineWidth',2);
    str = sprintf('Time: %d to %d ms from 0',time_win(i-2)-101,time_win(i-2)+50-101 );
    title(str);
    waitforbuttonpress;
end

figure;
time_win = [101:5:301];
for i = 3:43   
    temp = (a3.di_matrix{i,1} );% - a3.di_matrix{i,1});caxis([0 0.3]);
    temp(temp< thresh) =0;
    G = digraph(temp);
   plot(G,'XData',cell2mat(chanNamesLocs(:,2)),'YData',cell2mat(chanNamesLocs(:,3)),'LineWidth',2);%,'ZData',cell2mat(chanNamesLocs(:,4)));
    %plot(G,'layout','force','LineWidth',2);%,'ZData',cell2mat(chanNamesLocs(:,4)));
 
    str = sprintf('Time: %d to %d ms from 0',time_win(i-2)-101,time_win(i-2)+50-101 );
    title(str);
    waitforbuttonpress;
end
