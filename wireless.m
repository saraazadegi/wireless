clc;
clear;
close all;
%% part A  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% num_of_bits = 1000;
% bit_stream = randi([0,1],1,num_of_bits);
% symbol = 2*bit_stream - 1;
% Eb = 1;
% % SNRdB = 0;
% SNRdB = -15:5:35;
% SNR = 10.^(SNRdB/10);
% No = Eb ./ SNR;
% received = zeros(length(SNR),num_of_bits);
% symbol_error = zeros(1,length(SNR));
% for counter = 1:length(SNR)
%     N = sqrt(No(1,counter)/2)*randn(1,num_of_bits);
%     received(counter,:) = N + symbol;
%     scatterplot(received(counter,:));
%     title(strcat('SNR in dB is',num2str(SNRdB(1,counter))));
%     for i = 1:num_of_bits
%         if ((received(counter,i) > 0) && (bit_stream(1,i) == 0)) ||...
%                 ((received(counter,i) < 0) && (bit_stream(1,i) == 1))
%             symbol_error(1,counter) = symbol_error(1,counter) +1;
%         end
%     end
% end
% 
% symbol_error = symbol_error / num_of_bits;
% BER = symbol_error;
% figure;
% plot(SNRdB,BER,'--o');
% title('BER Vs SNR');
% xlabel('SNR per bir in dB');
% ylabel('BER%');
% scatterplot(symbol);
% title('sending symbols');


% %% part B  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eb = 1;
Es = 2;
SNRdB = -15:5:35; % Eb/No   SNR per bit
SNRdB = SNRdB * (Es/Eb); % Es/No SNR per Symbol
SNR = 10.^(SNRdB/10);
No = Es ./ SNR;
num_of_bits = 2000;
num_of_symbols = num_of_bits/2;
bit_stream = randi([0,1],1,num_of_bits);
symbol = zeros(1,num_of_symbols);
const1 = sqrt(Es)*[1i, 1, -1, -1i]; % bit mapping on symbols [00, 01, 10, 11]
const2 = sqrt(Es)*[1+1i, 1-1i, -1+1i, -1-1i]/sqrt(2); % [00, 01, 10, 11]

j = 0; % this variable using for switching between constelations
%% mappping bits on symbols
for i = 1:num_of_symbols
    if j  % pick symbol among const1
        if [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [0,0]
            symbol(1,i) = const1(1,1);
        elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [0,1]
            symbol(1,i) = const1(1,2);
        elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [1,0]
            symbol(1,i) = const1(1,3);
        elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [1,1]
            symbol(1,i) = const1(1,4);
        end
    else % pick symbol among const2 
        if [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [0,0]
            symbol(1,i) = const2(1,1);
        elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [0,1]
            symbol(1,i) = const2(1,2);
        elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [1,0]
            symbol(1,i) = const2(1,3);
        elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [1,1]
            symbol(1,i) = const2(1,4);
        end
    end
    j = ~j;
end

j = 0;
received = zeros(length(SNR),num_of_symbols);
received_bit_stream = zeros(length(SNR),num_of_bits);
BER = zeros(1,length(SNR));
for counter = 1:length(SNR)
    NI = sqrt(No(1,counter)/2)*randn(1,num_of_symbols);
    NQ = sqrt(No(1,counter)/2)*randn(1,num_of_symbols);
    received(counter,:) = NI + 1i*NQ + symbol;
    scatterplot(received(counter,:));
    title(strcat('SNR per bit in dB is',num2str(SNRdB(1,counter)/2)));
    for i = 1:num_of_symbols  %Detection
        if j  % pick symbol among const1
           if (angle(received(counter,i)) > pi/4) && ...
                   (angle(received(counter,i)) < pi*3/4)
                received_bit_stream(counter,2*i-1) = 0;
                received_bit_stream(counter,2*i) = 0;
            elseif (angle(received(counter,i)) > 3*pi/4) || ...
                    (angle(received(counter,i)) < -pi*3/4)
                received_bit_stream(counter,2*i-1) = 1;
                received_bit_stream(counter,2*i) = 0;
            elseif (angle(received(counter,i)) < pi/4) && ...
                    (angle(received(counter,i)) > -pi/4)
                received_bit_stream(counter,2*i-1) = 0;
                received_bit_stream(counter,2*i) = 1;
            elseif (angle(received(counter,i)) < -pi/4) && ...
                    (angle(received(counter,i)) > -pi*3/4)
                received_bit_stream(counter,2*i-1) = 1;
                received_bit_stream(counter,2*i) = 1;
            end
        else % pick symbol among const2
            if (real(received(counter,i)) > 0) && ...
                    (imag(received(counter,i)) > 0)
                received_bit_stream(counter,2*i-1) = 0;
                received_bit_stream(counter,2*i) = 0;
            elseif (real(received(counter,i)) < 0) && ...
                    (imag(received(counter,i)) > 0)
                received_bit_stream(counter,2*i-1) = 1;
                received_bit_stream(counter,2*i) = 0;
            elseif (real(received(counter,i)) > 0) && ...
                    (imag(received(counter,i)) < 0)
                received_bit_stream(counter,2*i-1) = 0;
                received_bit_stream(counter,2*i) = 1;
            elseif (real(received(counter,i)) < 0) && ...
                    (imag(received(counter,i)) < 0)
                received_bit_stream(counter,2*i-1) = 1;
                received_bit_stream(counter,2*i) = 1;
            end
        end
        j = ~j;
    end
end

for i = 1:length(SNR)
    BER(1,i) = sum(received_bit_stream(i,:) ~= bit_stream);
end
BER = BER/num_of_bits;
figure;
plot(SNRdB*Eb/Es,BER*100,'--o');
title('BER Vs SNR');
xlabel('SNR per bir in dB');
ylabel('BER%');

%% part C  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% part D  %%%%%%%%%%%%BPSK rayleigh channel effect%%%%%%%%%%%%%%%%%%%
% num_of_bits = 8000;
% bit_stream = randi([0,1],1,num_of_bits);
% symbol = 2*bit_stream - 1;
% Eb = 1;
% % SNRdB = 0;
% SNRdB = -15:5:35;
% SNR = 10.^(SNRdB/10);
% No = Eb ./ SNR;
% fading_rate = [0.0001, 0.001, 0.01, 0.1];
% fm = 10; % fm is the maximum Doppler shift
% % ts is the sample time of the input signal
% received = zeros(length(fading_rate),length(SNR),num_of_bits);
% symbol_error = zeros(length(fading_rate),length(SNR));
% k = 0;
% for ts = fading_rate/fm
%     k = k+1; % ts counter
%     ray_channel = rayleighchan(ts,fm);
%     symbol_after_ray = filter(ray_channel,symbol);
%     for counter = 1:length(SNR)
%         N = sqrt(No(1,counter)/2)*randn(1,num_of_bits);
%         received(k,counter,:) = N + symbol_after_ray;
% %         scatterplot(received(k,counter,:));
% %         title(strcat('SNR in dB is',num2str(SNRdB(1,counter)),...
% %             'fading rate is',num2str(fading_rate(1,k))));
%         for i = 1:num_of_bits
%             if ((received(k,counter,i) > 0) && (bit_stream(1,i) == 0)) ||...
%                     ((received(k,counter,i) < 0) && (bit_stream(1,i) == 1))
%                 symbol_error(k,counter) = symbol_error(k,counter) +1;
%             end
%         end
%     end
% end
% BER = symbol_error/num_of_bits;
% for k = 1:length(fading_rate)
%     figure;
%     plot(SNRdB,BER(k,:)*100,'--o');
%     title(strcat('BER Vs SNR for fading rate ',num2str(fading_rate(1,k))));
%     xlabel('SNR per bir in dB');
%     ylabel('BER%');
% end

%% part D  %%%%%%%%%%QPSK pi by 4 rayleigh channel effect%%%%%%%%%%%%%%%%
% 
% Eb = 1;
% Es = 2;
% SNRdB = -15:5:35; % Eb/No   SNR per bit
% SNRdB = SNRdB * (Es/Eb); % Es/No SNR per Symbol
% SNR = 10.^(SNRdB/10);
% No = Es ./ SNR;
% num_of_bits = 2000;
% num_of_symbols = num_of_bits/2;
% bit_stream = randi([0,1],1,num_of_bits);
% symbol = zeros(1,num_of_symbols);
% const1 = sqrt(Es)*[1i, 1, -1, -1i]; % bit mapping on symbols [00, 01, 10, 11]
% const2 = sqrt(Es)*[1+1i, 1-1i, -1+1i, -1-1i]/sqrt(2); % [00, 01, 10, 11]
% 
% j = 0; % this variable using for switching between constelations
% %% mappping bits on symbols
% for i = 1:num_of_symbols
%     if j  % pick symbol among const1
%         if [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [0,0]
%             symbol(1,i) = const1(1,1);
%         elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [0,1]
%             symbol(1,i) = const1(1,2);
%         elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [1,0]
%             symbol(1,i) = const1(1,3);
%         elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [1,1]
%             symbol(1,i) = const1(1,4);
%         end
%     else % pick symbol among const2 
%         if [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [0,0]
%             symbol(1,i) = const2(1,1);
%         elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [0,1]
%             symbol(1,i) = const2(1,2);
%         elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [1,0]
%             symbol(1,i) = const2(1,3);
%         elseif [bit_stream(1,2*i-1),bit_stream(1,2*i)] == [1,1]
%             symbol(1,i) = const2(1,4);
%         end
%     end
%     j = ~j;
% end
% j = 0;
% fading_rate = [0.0001, 0.001, 0.01, 0.1];
% fm = 10; % fm is the maximum Doppler shift
% % ts is the sample time of the input signal
% received_symbol = zeros(length(fading_rate),length(SNR),num_of_symbols);
% received_bit_stream = zeros(length(fading_rate),length(SNR),num_of_bits);
% symbol_error = zeros(length(fading_rate),length(SNR));
% BER = zeros(length(fading_rate),length(SNR));
% k = 0;
% for ts = fading_rate/fm
%     k = k+1; % ts counter
%     ray_channel = rayleighchan(ts,fm);
%     symbol_after_ray = filter(ray_channel,symbol);
%     for counter = 1:length(SNR)
%         NI = sqrt(No(1,counter)/2)*randn(1,num_of_symbols);
%         NQ = sqrt(No(1,counter)/2)*randn(1,num_of_symbols);
%         received_symbol(k,counter,:) = NI + 1i*NQ + symbol_after_ray;
%         scatterplot(received_symbol(k,counter,:));
%         title(strcat('SNR per bit in dB is',num2str(SNRdB(1,counter)*Eb/Es),...
%             'fading rate is',num2str(fading_rate(1,k))));
%         for i = 1:num_of_symbols  %Detection
%             if j  % pick symbol among const1
%                 if (angle(received_symbol(k,counter,i)) > pi/4) && ...
%                         (angle(received_symbol(k,counter,i)) < pi*3/4)
%                     received_bit_stream(k,counter,2*i-1) = 0;
%                     received_bit_stream(k,counter,2*i) = 0;
%                 elseif (angle(received_symbol(k,counter,i)) > 3*pi/4) || ...
%                         (angle(received_symbol(k,counter,i)) < -pi*3/4)
%                     received_bit_stream(k,counter,2*i-1) = 1;
%                     received_bit_stream(k,counter,2*i) = 0;
%                 elseif (angle(received_symbol(k,counter,i)) < pi/4) && ...
%                         (angle(received_symbol(k,counter,i)) > -pi/4)
%                     received_bit_stream(k,counter,2*i-1) = 0;
%                     received_bit_stream(k,counter,2*i) = 1;
%                 elseif (angle(received_symbol(k,counter,i)) < -pi/4) && ...
%                         (angle(received_symbol(k,counter,i)) > -pi*3/4)
%                     received_bit_stream(k,counter,2*i-1) = 1;
%                     received_bit_stream(k,counter,2*i) = 1;
%                 end
%             else % pick symbol among const2
%                 if (real(received_symbol(k,counter,i)) > 0) && ...
%                         (imag(received_symbol(k,counter,i)) > 0)
%                     received_bit_stream(k,counter,2*i-1) = 0;
%                     received_bit_stream(k,counter,2*i) = 0;
%                 elseif (real(received_symbol(k,counter,i)) < 0) && ...
%                         (imag(received_symbol(k,counter,i)) > 0)
%                     received_bit_stream(k,counter,2*i-1) = 1;
%                     received_bit_stream(k,counter,2*i) = 0;
%                 elseif (real(received_symbol(k,counter,i)) > 0) && ...
%                         (imag(received_symbol(k,counter,i)) < 0)
%                     received_bit_stream(k,counter,2*i-1) = 0;
%                     received_bit_stream(k,counter,2*i) = 1;
%                 elseif (real(received_symbol(k,counter,i)) < 0) && ...
%                         (imag(received_symbol(k,counter,i)) < 0)
%                     received_bit_stream(k,counter,2*i-1) = 1;
%                     received_bit_stream(k,counter,2*i) = 1;
%                 end
%             end
%             j = ~j;
%         end
%     end
% end
% for k = 1:length(fading_rate)
%     for i = 1:length(SNR)
%         BER(k,i) = sum(squeeze(received_bit_stream(k,i,:))' ~= bit_stream);
%     end
% end
% BER = BER/num_of_bits;
% for k = 1:length(fading_rate)
%     figure;
%     plot(SNRdB,BER(k,:)*100,'--o');
%     title(strcat('BER Vs SNR for fading rate ',num2str(fading_rate(1,k))));
%     xlabel('SNR per bir in dB');
%     ylabel('BER%');
% end