% symbol error rate for QPSK(4-QAM) modulation
out=0.1;
QPSK_x=zeros(3,0);
QPSK_y=zeros(3,0);
for i=1:3
    QPSK_F=qfuncinv(out)^2/2;
    QPSK_x(i)=10*log10(QPSK_F);
    QPSK_y(i)=out;
    out=out/10;
end
N = 10^5; % number of symbols
Es_N0_dB = QPSK_x; % multiple Eb/N0 values
ipHat = zeros(1,N);
for i = 1:length(Es_N0_dB)
    ip = (2*(rand(1,N)>0.5)-1) + 1i*(2*(rand(1,N)>0.5)-1); %
    s = (1/sqrt(2))*ip; % normalization of energy to 1
    n = 1/sqrt(2)*[randn(1,N) + 1i*randn(1,N)]; % white guassian noise, 0dB variance

    y = s + 10^(-Es_N0_dB(i)/20)*n; % additive white gaussian noise

    y_re = real(y); % real
    y_im = imag(y); % imaginary
    ipHat(find(y_re < 0 & y_im < 0)) = -1 + -1*1i;
    ipHat(find(y_re >= 0 & y_im > 0)) = 1 + 1*1i;
    ipHat(find(y_re < 0 & y_im >= 0)) = -1 + 1*1i;
    ipHat(find(y_re >= 0 & y_im < 0)) = 1 - 1*1i;
    nErr(i) = size(find([ip- ipHat]),2); % couting the number of errors
end

simSer_QPSK = nErr/N;

figure
semilogy(Es_N0_dB,QPSK_y,'b.-');
hold on
semilogy(Es_N0_dB,simSer_QPSK,'mx-');
grid on
legend('theory-QPSK', 'simulation-QPSK');
xlabel('Es/No, dB')
ylabel('Symbol Error Rate')
title('Symbol error probability curve for QPSK')


