clc
clear
close all
%% Data Acquisition 

%Import Data 
data = xlsread('Data');   %Fz|Cz|O2|HEOG

%Data Extraction
eog_HEOG = -data(:,4);     %HEOG artifact
eeg_Fz = data(:,1);       %Raw at Fz
eeg_Cz = data(:,2);       %Raw at Cz
eeg_O2 = data(:,3);       %Raw at O2

%Initial data plots
figure(1)

subplot(4,1,1)              %HEOG plot
plot(eog_HEOG)
title('Horizontal EOG Artifact')
xlabel('sample number')
ylabel('HEOG mu volts')

subplot(4,1,2)              %Raw at Fz plot
plot(eeg_Fz)
title('Raw EEG at Fz')
xlabel('sample number')
ylabel('EEG mu volts')

subplot(4,1,3)              %Raw at Cz plot
plot(eeg_Cz)
title('Raw EEG at Cz')
xlabel('sample number')
ylabel('EEG mu volts')

subplot(4,1,4)              %Raw at O2 plot
plot(eeg_O2)
title('Raw EEG at O2')
xlabel('sample number')
ylabel('EEG mu volts')
