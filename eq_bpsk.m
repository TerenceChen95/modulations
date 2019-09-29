N  = 10^6; % number of bits or symbols
Eb_N0_dB = [0:15]; % multiple Eb/N0 values
K = 4;

for i = 1:length(Eb_N0_dB)

   % Transmitter
   ip = rand(1,N)>0.5; % generating 0,1 with equal probability
   s = 2*ip-1; % BPSK modulation 0 -> -1; 1 -> 0 

   % Channel model, multipath channel
   nTap = 3;
   E=0.5^2+0.25^2+0.125^2;
   ht=[0.5/(sqrt(E)),0.25/sqrt(E),0.125/sqrt(E)];
   L  = length(ht);

   chanOut = conv(s,ht);  
   n = 1/sqrt(2)*[randn(1,N+length(ht)-1) + 1i*randn(1,N+length(ht)-1)]; % white gaussian noise, 0dB variance 
   
   % Noise addition
   y = chanOut + 10^(-Eb_N0_dB(i)/20)*n; % additive white gaussian noise 
   % mmse equalization
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
       %mmse
       hAutoCorr = conv(ht,fliplr(ht));
       hM = toeplitz([hAutoCorr([3:end]) zeros(1,2*k+1-L)], [ hAutoCorr([3:end]) zeros(1,2*k+1-L) ]);
       hM = hM + 1/2*10^(-Eb_N0_dB(i)/10)*eye(2*k+1);
       d  = zeros(1,2*k+1);
       d([-1:1]+k+1) = fliplr(ht);
       c_mmse  = [inv(hM)*d.'].';
       yFilt_mmse = conv(y,c_mmse);
       yFilt_mmse = yFilt_mmse(k+2:end); 
       yFilt_mmse = conv(yFilt_mmse,ones(1,1)); % convolution
       ySamp_mmse = yFilt_mmse(1:1:N);  % sampling at time T
       
       ipHat_mmse = real(ySamp_mmse)>0;
       nErr_mmse(k,i) = size(find([ip- ipHat_mmse]),2);
      % receiver - hard decision decoding
       ipHat_zf = real(ySamp_zf)>0;

       % counting the errors
       nErr_zf(k,i) = size(find([ip- ipHat_zf]),2);
   end   
end
simBer_zf = nErr_zf/N; % simulated ber
simBer_mmse = nErr_mmse/N; % simulated ber

%LMS  
%calculate delta
sum=0;
for i=1:N
    sum=sum+y(1,i);
end
expe_u=sum/N;
sum_2=0;
for i=1:N
    sum_2=sum_2+y(1,i)-expe_u;
end
delta=0.002/(expe_u^2+(sum_2/N)^2);
%
l=length(c_mmse);
t=[0:l-1];
d=s;
w=zeros(1,N);
e=zeros(1,N);
for i=1:N
    e(i)=d(i)-w(i)'*y(1,i);
    w(i+1)=w(i)+delta*e(i)*y(1,i);
end
dis=zeros(1,length(c_mmse));
for i=1:length(c_mmse)
    dis(i)=sqrt(abs(c_mmse(i)^2-w(i)^2));
end
figure,
plot(t,dis);
xlabel('t')
ylabel('Distance')
title('distance between MMSE and LMS')

%previous AWGN BER
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




% plot

figure
semilogy(Eb_N0_dB,simBer_zf(2,:),'bd-');
hold on
semilogy(Eb_N0_dB,simBer_zf(3,:),'bx-');
hold on
semilogy(Eb_N0_dB,simBer_zf(4,:),'bo-');
hold on
semilogy(Eb_N0_dB,simBer_mmse(2,:),'gd-');
hold on
semilogy(Eb_N0_dB,simBer_mmse(3,:),'gx-');
hold on
semilogy(Eb_N0_dB,simBer_mmse(4,:),'go-');
hold on 
semilogy(Eb_N0_dB,simBer,'mx-');
axis([0 14 10^-5 0.5])
grid on
legend('zf-5tap','zf-7tap','zf-9tap', 'mmse-5tap','mmse-7tap','mmse-9tap','AWGN');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for BPSK in ISI with ZF/MMSE equalizer');