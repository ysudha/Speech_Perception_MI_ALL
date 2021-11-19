function [data] = extract_VOT_data(subjCell, patient_number, VOT_bin)

% data is in subjCell size 21 x 2500 x 2
% extract all VOT data from patient number provided
trial =1; % loop through all trials
idx = 1; % index of trials for that VOT bin
data = [];

while ~isempty(subjCell{patient_number,trial,2}) % while the data is still there
    trial_id = subjCell{patient_number,trial,2}; % extract trial id
    if(trial_id(5) == VOT_bin)
        data{idx} = subjCell{patient_number,trial,1};
        idx = idx +1;
    end
    
    trial = trial +1;
end
        

end