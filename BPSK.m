out=0.1;
BPSK_x=zeros(3,1);
BPSK_y=zeros(3,1);
for i=1:3
    BPSK_F=(qfuncinv(out/0.5)^2)/2;
    BPSK_x(i)=10*log10(BPSK_F);
    BPSK_y(i)=out;
    out=out/10;
end
Eb_N0_dB=BPSK_x;
N = 10^5; % number of symbols
ipHat = zeros(1,N);
for i=1:length(Eb_N0_dB)
    input=rand(1,N)>0.5;   %generate random inputs
    s=2*input-1;    %BPSK modulation
    n = 1/sqrt(2)*[randn(1,N) + 1i*randn(1,N)]; % white guassian noise, 0dB variance
    y = s + 10^(-Eb_N0_dB(i)/20)*n; % additive white gaussian noise

    ipHat = real(y)>0;
    nErr(i) = size(find([input- ipHat]),2); % couting the number of errors   
end
simBer = nErr/N;
figure
%semilogy(Eb_N0_dB,BPSK_y,'b.-');
%hold on
semilogy(Eb_N0_dB,simBer,'mx-');
grid on
legend('simulation');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for BPSK modulation')
