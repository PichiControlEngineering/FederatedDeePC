fs_ctrl = 250; % ctrl sample rate
inputSignal_directory = convertCharsToStrings(findFolder("inputSignals"));
multisine_dir = inputSignal_directory+"\MultiSin2_"+fs_ctrl+"Hz.mat";

% % GENERATE DATA (optional)
% fmin = 5; fmax = fs_ctrl/2;
% sig = MultiSineGenerator(1+fs_ctrl*24, fmin, fmax, fs_ctrl);
% inputData = timetable(sig/max(sig)*2, 'SampleRate', 250);
% save(multisine_dir, "inputData")

%% 
load(multisine_dir)