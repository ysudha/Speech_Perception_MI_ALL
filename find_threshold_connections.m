function [thresh p_thresh] = find_threshold_connections(all_DI_values, art_win_num, top_per)
nch = size(all_DI_values,2);
a = all_DI_values(art_win_num,:,:);
b = a(:);
b = sort(b,'descend');
n = ceil(top_per*length(b)/100); % top x percent of connections threshold
thresh = b(n);

p_thresh=[];
for i = 1:size(all_DI_values,1)
    a = all_DI_values(i,:,:);
    
    a(a<thresh)=0;
    a = squeeze(a);
    p_thresh(i) = nnz(a(:));
end

p_thresh = p_thresh/(nch*nch-nch);
end