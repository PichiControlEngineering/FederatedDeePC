%% After running experiments, use this to save
PatoNR = 10; % Current Pato Number

%
ExpData_directory = convertCharsToStrings(findFolder("Experiment Data"));
data = logsout;
save_file_directory = ExpData_directory+"\"+ PatoNR+"\PatoResp2_"+fs_ctrl+"Hz.mat";

% Ensure target folder exists before saving
save_dir = fileparts(save_file_directory); % get directory portion
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

save(save_file_directory, "data")
disp("Data saved succesfully")