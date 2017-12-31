clc
clear
close all
%% Data Acquisition 

%Import Data 
data = xlsread('Data');   %Fz|Cz|Oz|HEOG

%Data Extraction
eog_HEOG = -data(:,4);     %HEOG artifact
eeg_Fz = data(:,1);       %Raw at Fz
eeg_Cz = data(:,2);       %Raw at Cz
eeg_O2 = data(:,3);       %Raw at O2

%% RLS ALGORITHM EXECUTION FOR Denoising at Fz 

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

%% RLS ALGORITHM EXECUTION FOR Denoising at Cz 
 

% Execution for Cz electrode
% Variable Initialization
sample_no = size(eeg_Cz);    % No of samples/time points
order = 3;                   % Order of the Adaptive Filter (User Tunable)
sigma = 0.01;                % Initializing variable
lambda = 0.9999;             % Forgetting Factor for RLS Algorithm (User Tunable)
H = zeros(order,1);          % Initial filter Coefficients
R = sigma*eye(order,order);  % Initial value for Reference combination


% RLS Algorithm for Adaptive Denoising
for n = 1:sample_no          % Loop to simulate reality situation
   
    s = eeg_Cz(n,1);         % eeg @ Cz at that time point
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
        correct_Cz(n,1) = s - (r'*H);
        
    end
end

% Cz Denoising plot
figure(3)

subplot(3,1,1)              %HEOG plot
plot(eog_HEOG)
title('Horizontal EOG Artifact')
xlabel('sample number')
ylabel('HEOG mu volts')

subplot(3,1,2)              %Raw at Cz plot
plot(eeg_Cz)
title('Raw EEG at Cz')
xlabel('sample number')
ylabel('EEG mu volts')

subplot(3,1,3)              %denoised eeg at Cz plot
plot(correct_Cz)
title('Corrected EEG at Cz')
xlabel('sample number')
ylabel('EEG mu volts')


%% RLS ALGORITHM EXECUTION FOR Denoising at O2 

% Execution for O2 electrode
% Variable Initialization
sample_no = size(eeg_O2);    % No of samples/time points
order = 3;                   % Order of the Adaptive Filter (User Tunable)
sigma = 0.01;                % Initializing variable
lambda = 0.9999;             % Forgetting Factor for RLS Algorithm (User Tunable)
H = zeros(order,1);          % Initial filter Coefficients
R = sigma*eye(order,order);  % Initial value for Reference combination


% RLS Algorithm for Adaptive Denoising
for n = 1:sample_no          % Loop to simulate reality situation
   
    s = eeg_O2(n,1);         % eeg @ O2 at that time point
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
        correct_O2(n,1) = s - (r'*H);
        
    end
end

% Oz Denoising plot
figure(4)

subplot(3,1,1)              %HEOG plot
plot(eog_HEOG)
title('Horizontal Noise Artifact')
xlabel('sample number')
ylabel('HEOG mu volts')

subplot(3,1,2)              %Raw at O2 plot
plot(eeg_O2)
title('Raw EEG at O2')
xlabel('sample number')
ylabel('EEG mu volts')

subplot(3,1,3)              %denoised eeg at O2 plot
plot(correct_O2)
title('Corrected EEG at O2')
xlabel('sample number')
ylabel('EEG mu volts')

