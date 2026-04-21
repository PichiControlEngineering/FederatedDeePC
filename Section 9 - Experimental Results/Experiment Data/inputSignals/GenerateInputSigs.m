%% Generate Input signals for the PATO
path = "C:\Users\20192556\Documents\Y6\4. Experiments\SpeedGoatExperiment\inputSignals";
% We generate a Multisin at 100, 250 and 1000 Hz, stored at 
fs = 1e3; Ts = 1./fs;

% Setting freq range for multisin
f_min = 2.5; f_max = 150;

fs_des = [100, 250, 1000];
T_test = 120; % length of data-measurement is seconds

for fs_i = fs_des
% Generate multisine input signal
u = MultiSineGenerator(T_test*fs_i,f_min,f_max,fs_i);
u = 2.4.*u./max(u); % normalizing
pspectrum(u,fs)

% inputData = array2timetable(u,"SampleRate",fs);
inputData = array2timetable(u,"SampleRate",fs_i); %downsampled data

% only add new input signals if they don't exist yet
files = dir(path);
files = files(~ismember({files.name}, {'.', '..'})); % Remove '.' and '..' entries

if isempty(files)
    save(strcat(path,"\MultiSin_",num2str(fs_i),"Hz"), "inputData", "fs_i")
end
end