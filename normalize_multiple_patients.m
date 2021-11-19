function bb = normalize_multiple_patients(aa)
% aa is a 21 x 3 matrix, 21 patients, 3 conditions
global_mean = mean(aa(:));
mean_each_row = mean(aa,2); % mean across 2nd dimension, for each row

bb = aa - repmat(mean_each_row,1,size(aa,2));
bb = bb + global_mean;

end