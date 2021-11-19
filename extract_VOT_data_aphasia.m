function [data] = extract_VOT_data_aphasia(subjCell, VOT_bin)

% data is in subjCell size 30 x 360 x 1377
% extract all VOT data from patient number provided
idx = 1; % index of trials for that VOT bin
data = [];

for trial = 1: size(subjCell,3) % while the data is still there
    trial_id = subjCell(1,5,trial); % extract trial id
    if(trial_id == VOT_bin)
        data{idx} = subjCell(:,11:end,trial);
        idx = idx +1;
    end

end
        

end