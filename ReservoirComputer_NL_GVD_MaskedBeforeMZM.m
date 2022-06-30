clc;
close all;
clear;

[Input,Output] = loadDataset();
[Pulse,Fiber] = loadSimulationParameters();

%% test split step function
[SSFM_t,SSFM_f] = splitStepMethod(Pulse.guassian,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega );
plot(Pulse.Normalize_tau,abs(Pulse.guassian).^2);hold on;
plot(Pulse.Normalize_tau,abs(SSFM_t).^2);

%% Reservoir Network application to NARMA10 task 
input_mask = [0.7 0.2 0.4 0.5 0.9 0.8 0.1	0.4	0.3	0.6	0.95	0.3	0.99	0.35	0.2	0.7	0.8	0.6	0.1	0.4	0.2	0.2	0.1	0.9	0.8	0.4	0.8	0.75	0.6	0.4	0.35	0.15	0.3	0.55	0.45	0.65	0.5	0.8	0.1	0.4	0.6	0.2	0.6	0.1	0.8	0.7	0.3	0.4	1.0	0.3];
u_n = Input(1:1000);
Y_n = Output(1:1000);

mask_len = length(input_mask);
signal_len = length(u_n);
    
reservoirmemory = zeros(mask_len,Pulse.fft_point);
readout_t = zeros(signal_len,Pulse.fft_point);
readout_f = zeros(signal_len,Pulse.fft_point);
readout_RS = zeros(signal_len, mask_len);
%% Algorithm test, Single delay loop with T = (N+1)Theta
% for i = 1:signal_len
%     time_multiplexed_signal = input_mask .* u_n(i);    
%     for k = 1:length(time_multiplexed_signal)-1   %T = (N+1)Theta
%         modulated_signal = time_multiplexed_signal(k) .* Pulse.guassian;
%         memory_modulated_signal = modulated_signal + reservoirmemory(k+1,:);
%         [SSFM_t,SSFM_f] = splitStepMethod(memory_modulated_signal,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega ); %split step algorithm to resolve GVD and NL effect in the fiber
%         SSFM_f =  SSFM_f .* (sqrt(2*pi)/(Pulse.fft_point*Pulse.delta_tau));
%         reservoirmemory(k+1,:) = SSFM_t;        
%         %figure; plot(Pulse.Normalize_tau,abs(reservoirmemory(k,:)).^2);
%     end
%     %%last data
%     modulated_signal = time_multiplexed_signal(length(time_multiplexed_signal)) .* Pulse.guassian;
%     memory_modulated_signal = modulated_signal + reservoirmemory(1,:);
%     [SSFM_t,SSFM_f] = splitStepMethod(memory_modulated_signal,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega ); %split step algorithm to resolve GVD and NL effect in the fiber
%     SSFM_f =  SSFM_f .* (sqrt(2*pi)/(Pulse.fft_point*Pulse.delta_tau));
%     reservoirmemory(1,:) = SSFM_t;
%     readout_t(i,:) = abs(reservoirmemory(round(k/2),:)).^2;
%     %figure; plot(Pulse.Normalize_tau,readout_t(i,:));
% end  

%% Multiple Delay algorithm according to paper
% input_mask2 = [0.910000000000000	0.130000000000000	0.910000000000000	0.630000000000000	0.100000000000000	0.280000000000000	0.550000000000000	0.960000000000000	0.960000000000000	0.160000000000000	0.970000000000000	0.960000000000000	0.490000000000000	0.800000000000000	0.140000000000000	0.420000000000000	0.920000000000000	0.790000000000000	0.960000000000000	0.660000000000000	0.0400000000000000	0.850000000000000	0.930000000000000	0.680000000000000	0.760000000000000	0.740000000000000	0.390000000000000	0.660000000000000	0.170000000000000	0.710000000000000	0.0300000000000000	0.280000000000000	0.0500000000000000	0.100000000000000	0.820000000000000	0.690000000000000	0.320000000000000	0.950000000000000	0.0300000000000000	0.440000000000000	0.380000000000000	0.770000000000000	0.800000000000000	0.190000000000000	0.490000000000000	0.450000000000000	0.650000000000000	0.710000000000000	0.750000000000000	0.280000000000000];
% delayVal = 5;
% %delayed_u_n = [zeros(1,delayVal),u_n(1:signal_len-delayVal)];
% % Algorithm test
% for i = delayVal+1:signal_len
%     time_multiplexed_signal = (input_mask .* u_n(i)) + (input_mask2 .* u_n(i-delayVal));    
%     for k = 1:length(time_multiplexed_signal)-1   %T = (N+1)Theta
%         modulated_signal = time_multiplexed_signal(k) .* Pulse.guassian;
%         memory_modulated_signal = modulated_signal + reservoirmemory(k+1,:);
%         [SSFM_t,SSFM_f] = splitStepMethod(memory_modulated_signal,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega ); %split step algorithm to resolve GVD and NL effect in the fiber
%         %SSFM_f =  SSFM_f .* (sqrt(2*pi)/(Pulse.fft_point*Pulse.delta_tau));
%         reservoirmemory(k+1,:) = SSFM_t;        
%         %figure; plot(Pulse.Normalize_tau,abs(reservoirmemory(k,:)).^2);
%     end
%     %%last data
%     modulated_signal = time_multiplexed_signal(length(time_multiplexed_signal)) .* Pulse.guassian;
%     memory_modulated_signal = modulated_signal + reservoirmemory(1,:);
%     [SSFM_t,SSFM_f] = splitStepMethod(memory_modulated_signal,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega ); %split step algorithm to resolve GVD and NL effect in the fiber
%     SSFM_f =  SSFM_f .* (sqrt(2*pi)/(Pulse.fft_point*Pulse.delta_tau));
%     reservoirmemory(1,:) = SSFM_t;
%     signalIntensity = abs(reservoirmemory(round(k/2),:)).^2;
%     readout_t(i,:) = signalIntensity ./max(signalIntensity);
%     %figure; plot(Pulse.Normalize_tau,readout_t(i,:));
% end
%% Algorithm for adjusting Reservoir state values
for i = 1:signal_len
    time_multiplexed_signal = input_mask .* u_n(i);    
    for k = 1:length(time_multiplexed_signal)-1   %T = (N+1)Theta
        modulated_signal = time_multiplexed_signal(k) .* Pulse.guassian;
        memory_modulated_signal = modulated_signal + reservoirmemory(k+1,:);
        [SSFM_t,SSFM_f] = splitStepMethod(memory_modulated_signal,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega ); %split step algorithm to resolve GVD and NL effect in the fiber
        %SSFM_f =  SSFM_f .* (sqrt(2*pi)/(Pulse.fft_point*Pulse.delta_tau));
        reservoirmemory(k+1,:) = SSFM_t;        
        readout_RS(i,k) = max(abs(SSFM_t).^2);
    end
    %%last data
    modulated_signal = time_multiplexed_signal(length(time_multiplexed_signal)) .* Pulse.guassian;
    memory_modulated_signal = modulated_signal + reservoirmemory(1,:);
    [SSFM_t,SSFM_f] = splitStepMethod(memory_modulated_signal,Pulse.soliton_order,Fiber.length, Fiber.Beta2, Fiber.Division_factor, Pulse.omega ); %split step algorithm to resolve GVD and NL effect in the fiber
    %SSFM_f =  SSFM_f .* (sqrt(2*pi)/(Pulse.fft_point*Pulse.delta_tau));
    reservoirmemory(1,:) = SSFM_t;
    readout_RS(i,50) = max(abs(SSFM_t).^2);
end  


%% output data to file
writematrix(readout_RS,'NARMA_TASKIN_RS.csv');
writematrix(Y_n,'NARMA_TASKOUT_RS.csv');