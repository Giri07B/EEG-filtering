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
eeg_Oz = data(:,3);       %Raw at Oz

%% Claim 5: As order increases the MSE monotonically decreases

% Execution for Fz electrode with M = 3
% Variable Initialization
sample_no = size(eeg_Fz);    % No of samples/time points
order = 3;                   % Order of the Adaptive Filter (User Tunable)
sigma = 0.01;                % Initializing variable
lambda = 0.9999;             % Forgetting Factor for RLS Algorithm (User Tunable)
H = zeros(order,1);          % Initial filter Coefficients
R = sigma*eye(order,order);  % Initial value for Reference combination

MSE_M1 = 0;
MSE_M2 = 0;
MSE_M3 = 0;


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
        correct_Fz_M1(n,1) = s - (r'*H);
        
        MSE_M1 = ((lambda^(sample_no(1,1)-n))*((correct_Fz_M1(n,1))^2))+ MSE_M1;
    end
end

MSE_M1 = MSE_M1/(sample_no(1,1) - order + 1)

% Execution for Fz electrode with M = 1
% Variable Initialization
sample_no = size(eeg_Fz);    % No of samples/time points
order = 1;                   % Order of the Adaptive Filter (User Tunable)
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
        
            r1(j,1) = eog_HEOG(i,1);
            j = j+1;
        end
        
        % Calculation of Filter Coefficients
        K = ((inv(R))*r1)*(inv(lambda+r1'*(inv(R))*r1)); 
        e = s - (r1'*H);
        H = H + (K*e);
        R = inv((inv(lambda)*inv(R)) - ((inv(lambda))*(K)*((r1)')*(inv(R))));
        % Calculation Of RLS algorithm estimated signal
        correct_Fz_M2(n,1) = s - (r1'*H);
     
        MSE_M2 = ((lambda^(sample_no(1,1)-n))*((correct_Fz_M2(n,1))^2))+ MSE_M2;
        
    end
end

MSE_M2 = MSE_M2/(sample_no(1,1) - order + 1)

% Execution for Fz electrode with M = 12
% Variable Initialization
sample_no = size(eeg_Fz);    % No of samples/time points
order = 12;                   % Order of the Adaptive Filter (User Tunable)
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
        correct_Fz_M3(n,1) = s - (r'*H);
         
        MSE_M3 = ((lambda^(sample_no(1,1)-n))*((correct_Fz_M3(n,1))^2))+ MSE_M3;
    
    end
end

MSE_M3= MSE_M3/(sample_no(1,1) - order + 1)
