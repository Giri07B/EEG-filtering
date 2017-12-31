clc
clear
close all
%% Data Acquisition 

%Import Data 
data = xlsread('Data');   %Fz|Cz|O2|HEOG

%Data Extraction
eog_HEOG = -data(:,4);     %HEOG artifact
eeg_Fz = data(:,1);       %Raw at Fz

%% RLS ALGORITHM EXECUTION FOR Denoising at Fz is a easy, stable and fast convergence method suitable for online removal of EOG artifacts 

% Execution for Fz electrode
% Variable Initialization
sample_no = size(eeg_Fz);    % No of samples/time points
order = 3;                   % Order of the Adaptive Filter (User Tunable)
sigma = 0.01;                % Initializing variable
lambda = 0.9999;             % Forgetting Factor for RLS Algorithm (User Tunable)
H = zeros(order,1);          % Initial filter Coefficients
R = sigma*eye(order,order);  % Initial value for Reference combination


% RLS Algorithm for Adaptive Denoising
for n = 1:sample_no          % Loop to simulate reality situation
   
    s = eeg_Fz(n,1);         % eeg @ Fz at that time point
    if n>=order               
            
        j = 1;               % calculation of last "order" reference vector
        for i = n:-1:(n+1-order)  
        
            r(j,1) = eog_HEOG(i,1);
            j = j+1;
        end
        
        % Calculation of Filter Coefficients
        K = ((inv(R))*r)*(inv(lambda+r'*(inv(R))*r)); 
        e = s - (r'*H);
        H = H + (K*e);
        R = inv((inv(lambda)*inv(R)) - ((inv(lambda))*(K)*((r)')*(inv(R))));
        % Calculation Of RLS algorithm estimated signal
        correct_Fz(n,1) = s - (r'*H);
        
    end
end

% Fz Denoising plot
figure(2)

subplot(3,1,1)              %HEOG plot
plot(eog_HEOG)
title('Horizontal EOG Artifact')
xlabel('sample number')
ylabel('HEOG mu volts')

subplot(3,1,2)              %Raw at Fz plot
plot(eeg_Fz)
title('Raw EEG at Fz')
xlabel('sample number')
ylabel('EEG mu volts')

subplot(3,1,3)              %denoised eeg at Fz plot
plot(correct_Fz)
title('Corrected EEG at Fz')
xlabel('sample number')
ylabel('EEG mu volts')

