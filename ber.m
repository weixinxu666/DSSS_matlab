clc;
clear all;
%% parameters
Fs = 1000;
fc = 100;
fp = 2;
bit_t = 0.1;
wq = 221;
nsnr = 25;
navg = 100;
snr_step = -2;
bit_error_1 = zeros(3,nsnr);
bit_error_2 = zeros(3,nsnr);
bit_error_3 = zeros(3,nsnr);
chan2 = ricianchan(1/Fs,0,3);
chan2.PathDelays = [0 1e-6 1e-7];
chan2.ResetBeforeFiltering = 0;
chan2.NormalizePathGains = 1;

for w = 1:3
%% message generation with BPSK
m = randi([0 1], 1, 20);
for bit = 1:length(m)
    if(m(bit)==0)
        m(bit) = -1;
    end
end
message = repmat(m,fp,1);
message = reshape(message,1,[]);
%% PN generation and multiply with message
pn_code = randi([0,1],1,length(m)*fp);
for bit = 1:length(pn_code)
    if(pn_code(bit)==0)
        pn_code(bit) = -1;
    end
end

DSSS = message.*pn_code;
%% create carrier and multipy with encoded sequence
t = 0:1/Fs:(bit_t-1/Fs);
s0 = -1*cos(2*pi*fc*t);
s1 = cos(2*pi*fc*t);
carrier = [];
BPSK = [];
for i = 1:length(DSSS)
    if (DSSS(i) == 1)
        BPSK = [BPSK s1];
    elseif (DSSS(i) == -1)
        BPSK = [BPSK s0];
    end
    carrier = [carrier s1];
end

 

rice = filter(chan2,BPSK);

for h=1:nsnr
    noise = [];
    BER_1 = 0;
    BER_2 = 0;
    BER_3 = 0;

    for y=1:navg
        noise = awgn(BPSK,(snr_step*h));
        rice_rx = awgn(rice,(snr_step*h));

        % Noise variance
        N0 = 1/10^((snr_step*h+11)/10);
        % Rayleigh channel fading
        ray_chan = 1/sqrt(2)*(randn(1,length(BPSK)) +1j*randn(1,length(BPSK)));
%         ray_chan = 1/sqrt(2)*randn(1,length(BPSK));
        % Send over Gaussian Link to the receiver
        rayl = ray_chan.*BPSK + sqrt(N0/2)*(randn(1,length(BPSK))+i*randn(1,length(BPSK)));
        % Equalization to remove fading effects. Ideal Equalization Considered
        rayl_rx = rayl./ray_chan;
        rx_1 =[];
        rx_2 =[];
        rx_3 =[];
        for i = 1:length(pn_code)
            if(pn_code(i)==1)
                rx_1 = [rx_1 noise((((i-1)*length(t))+1):i*length(t))];
                rx_2 = [rx_2 rayl_rx((((i-1)*length(t))+1):i*length(t))];
                rx_3 = [rx_3 rice_rx((((i-1)*length(t))+1):i*length(t))];
            else
                rx_1 = [rx_1 (-1)*noise((((i-1)*length(t))+1):i*length(t))];
                rx_2 = [rx_2 (-1)*rayl_rx((((i-1)*length(t))+1):i*length(t))];
                rx_3 = [rx_3 (-1)*rice_rx((((i-1)*length(t))+1):i*length(t))];
            end
        end
        result_1 = zeros(1,length(m));
        result_2 = zeros(1,length(m));
        result_3 = zeros(1,length(m));
        for i = 1:length(m)
            x = length(t)*fp;
            cx_1 = sum(carrier(((i-1)*x)+1:i*x).*rx_1(((i-1)*x)+1:i*x));
            cx_2 = sum(carrier(((i-1)*x)+1:i*x).*rx_2(((i-1)*x)+1:i*x));
            cx_3 = sum(carrier(((i-1)*x)+1:i*x).*rx_3(((i-1)*x)+1:i*x));
            if(cx_1>0)
                result_1(i) = 1;
            else
                result_1(i) = -1;
            end
            if(cx_2>0)
                result_2(i) = 1;
            else
                 result_2(i) = -1;
            end           
            if(cx_3>0)
                result_3(i) = 1;
            else
                result_3(i) = -1;
            end
        end
        counter_1 = 0;
        counter_2 = 0;
        counter_3 = 0;
        for z=1:length(m)
            if m(z)~=result_1(z)
                counter_1 = counter_1+1;
            end

            if m(z)~=result_2(z)
                counter_2 = counter_2+1;
            end
            if m(z)~=result_3(z)
                counter_3 = counter_3+1;
            end
        end
        BER_1 = BER_1 + counter_1/length(m);
        BER_2 = BER_2 + counter_2/length(m);
        BER_3 = BER_3 + counter_3/length(m);
    end
    bit_error_1(w,h) = (BER_1/navg);
    bit_error_2(w,h) = (BER_2/navg);
    bit_error_3(w,h) = (BER_3/navg);
end
fp = fp+2;
end
 

%% draw BER - SNR graph
time_vec = (10):snr_step:((nsnr-1)*snr_step+10);
figure;
subplot(221)
semilogy(time_vec,bit_error_1(1,:),'g--o',time_vec,bit_error_1(2,:),'r--o',time_vec,bit_error_1(3,:),'b--o')
legend('2bits', '4bits','6bits')
axis([ -40 10 0 1 ])
title('BER for awgn noise');
xlabel('SNR values in dB');
ylabel('BER ratio');
subplot(222)
semilogy(time_vec,bit_error_2(1,:),'g--o',time_vec,bit_error_2(2,:),'r--o',time_vec,bit_error_2(3,:),'b--o')
legend('2bits', '4bits','6bits')
axis([ -30 30 0 1 ])
title('BER for Rayleigh channel');
xlabel('SNR values in dB');
ylabel('BER ratio');
subplot(223)
semilogy(time_vec,bit_error_3(1,:),'g--o',time_vec,bit_error_3(2,:),'r--o',time_vec,bit_error_3(3,:),'b--o')
legend('2bits', '4bits','6bits')
axis([ -40 10 0 1 ])
title('BER for Rician channel');
xlabel('SNR values in dB');
ylabel('BER ratio');
subplot(224)
semilogy(time_vec,bit_error_1(3,:),'g--o',time_vec,bit_error_2(3,:),'r--o',time_vec,bit_error_3(3,:),'b--o')
legend('BER of AWGN', 'BER of Rayleigh','BER for Rician')
axis([ -30 30 0 1 ])
title('BER for 6bits pn code ');
xlabel('SNR values in dB');
ylabel('BER ratio');
