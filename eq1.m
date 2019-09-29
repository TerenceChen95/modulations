%zero forcing for BPSK
N=10^6;
Eb_N0_dB=6.17200395505916; %SNR for BPSK to achieve 10^(-3) BER
K=5;
%Transmitter
input=rand(1,N)>0.5;   %generate random inputs
s=2*input-1;    %BPSK modulation
%channel model
mTap=5;
E=0.5^2+0.25^2+0.125^2;
ht=[0.5/(sqrt(E)),0.25/sqrt(E),0.125/sqrt(E)];
L  = length(ht);
chanOut=conv(s,ht);
n=1/sqrt(2)*[randn(1,N+length(ht)-1)+1i*randn(1,N+length(ht)-1)]; %white guassian noise
%add noise to the channel
y=chanOut+10^(-Eb_N0_dB/20)*n;
for k=2:K
% zero forcing equalization
hM = toeplitz([ht([2:end]) zeros(1,2*k+1-L+1)], [ ht([2:-1:1]) zeros(1,2*k+1-L+1) ]);
d  = zeros(1,2*k+1);
d(k+1) = 1;
c_zf  = [inv(hM)*d.'].';
yFilt_zf = conv(y,c_zf);
yFilt_zf = yFilt_zf(k+2:end); 
yFilt_zf = conv(yFilt_zf,ones(1,1)); % convolution
ySamp_zf = yFilt_zf(1:1:N);  % sampling at time T

% mmse equalization
hAutoCorr = conv(ht,fliplr(ht));
hM = toeplitz([hAutoCorr([3:end]) zeros(1,2*k+1-L)], [ hAutoCorr([3:end]) zeros(1,2*k+1-L) ]);
hM = hM + 1/2*10^(-Eb_N0_dB/10)*eye(2*k+1);
d  = zeros(1,2*k+1);
d([-1:1]+k+1) = fliplr(ht);
c_mmse  = [inv(hM)*d.'].';
yFilt_mmse = conv(y,c_mmse);
yFilt_mmse = yFilt_mmse(k+2:end); 
yFilt_mmse = conv(yFilt_mmse,ones(1,1)); % convolution
ySamp_mmse = yFilt_mmse(1:1:N);  % sampling at time T

% receiver - hard decision decoding
ipHat_zf = real(ySamp_zf)>0;
ipHat_mmse = real(ySamp_mmse)>0;

% counting the errors
nErr_zf(k) = size(find([input- ipHat_zf]),2);
nErr_mmse(k) = size(find([input- ipHat_mmse]),2);
end
simBer_zf = nErr_zf/N; % simulated ber
simBer_mmse = nErr_mmse/N; % simulated ber
theoryBer = 0.5*erfc(sqrt(10.^(Eb_N0_dB/10))); % theoretical ber

% plot
figure
semilogy(Eb_N0_dB,simBer_zf(:,1),'bs-');
hold on
semilogy(Eb_N0_dB,simBer_zf(:,2),'bx-');
semilogy(Eb_N0_dB,simBer_zf(:,3),'bd-');
semilogy(Eb_N0_dB,simBer_zf(:,4),'bo-');
semilogy(Eb_N0_dB,simBer_mmse(:,1),'gs-');
semilogy(Eb_N0_dB,simBer_mmse(:,2),'gx-');
semilogy(Eb_N0_dB,simBer_mmse(:,3),'gd-');
semilogy(Eb_N0_dB,simBer_mmse(:,4),'go-');
axis([0 14 10^-5 1])
grid on
legend('zf-5tap', 'zf-7tap','zf-9tap','zf-11tap','mmse-5tap','mmse-7tap','mmse-9tap','mmse-11tap')
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK in ISI with ZF/MMSE equalizer');


