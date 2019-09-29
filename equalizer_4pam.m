%zero forcing for BPSK
N=10^6;
Eb_N0_dB=[0:15]; %SNR for 4-PAM to achieve 10^(-3) BER
alpha4PAM=[-3,-1,1,3];
K=4;
for i = 1:length(Eb_N0_dB)
    %Transmitter
    input=randsrc(1,N,alpha4PAM);   %generate random inputs
    s=1/sqrt(5)*input;    %4PAM modulation

    %channel model
    mTap=5;
    E=0.5^2+0.25^2+0.125^2;
    ht=[0.5/(sqrt(E)),0.25/sqrt(E),0.125/sqrt(E)];  %normalize channel
    chanOut=conv(s,ht);
    n=1/sqrt(2)*[randn(1,N+length(ht)-1)+1i*randn(1,N+length(ht)-1)]; %white guassian noise

    %add noise to the channel
    y=chanOut+10^(-Eb_N0_dB(i)/20)*n;   
    L=length(ht);    
    for k=2:K
        %zero forcing
        hM=toeplitz([ht([2:end]) zeros(1,2*k+1-L+1)], [ht([2:-1:1]) zeros(1,2*k+1-L+1)]);
        d=zeros(1,2*k+1);
        d(k+1)=1;
        c=[inv(hM)*d.'].';
        %matched filter
        yF=conv(y,c);
        yF=yF(k+2:end);
        yF=conv(yF,ones(1,1));
        ySamp=yF(1:1:N);
        
        %mmse
        hAutoCorr = conv(ht,fliplr(ht));
        hM = toeplitz([hAutoCorr([3:end]) zeros(1,2*k+1-L)], [ hAutoCorr([3:end]) zeros(1,2*k+1-L) ]);
        hM = hM + 1/2*10^(-Eb_N0_dB(i)/10)*eye(2*k+1);
        d  = zeros(1,2*k+1);
        d([-1:1]+k+1) = fliplr(ht);
        c_mmse  = [inv(hM)*d.'].';
        yF_mmse = conv(y,c_mmse);
        yF_mmse = yF_mmse(k+2:end); 
        yF_mmse = conv(yF_mmse,ones(1,1)); % convolution
        ySamp_mmse = yF_mmse(1:1:N);  % sampling at time T

        %reciever
        r_mmse=real(ySamp_mmse);
        inputHat_mmse(find(r_mmse< -2/sqrt(5))) = -3;
        inputHat_mmse(find(r_mmse>= 2/sqrt(5))) = 3;
        inputHat_mmse(find(r_mmse>=-2/sqrt(5) & r_mmse<0)) = -1;
        inputHat_mmse(find(r_mmse>=0 & r_mmse<2/sqrt(5))) = 1;
        nErr_mmse(k,i) = size(find([input- inputHat_mmse]),2);
    end
    r_zf=real(ySamp);
    inputHat(find(r_zf< -2/sqrt(5))) = -3;
    inputHat(find(r_zf>= 2/sqrt(5))) = 3;
    inputHat(find(r_zf>=-2/sqrt(5) & r_zf<0)) = -1;
    inputHat(find(r_zf>=0 & r_zf<2/sqrt(5))) = 1;
    nErr_zf(k,i)=size(find([input-inputHat]),2);
end
simBer_zf = nErr_zf/N; % simulated ber ZF
simBer_mmse = nErr_mmse/N; % simulated ber MMSE

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

%previous AWGN
for i = 1:length(Eb_N0_dB)
    ip = randsrc(1,N,alpha4PAM);
    s = (1/sqrt(5))*ip; % normalization of energy to 1
    n = 1/sqrt(2)*[randn(1,N) + 1i*randn(1,N)]; % white guassian noise, 0dB variance
    y = s + 10^(-Eb_N0_dB(i)/20)*n; % additive white gaussian noise
    % demodulation
    r = real(y); % taking only the real part
    ipHat(find(r< -2/sqrt(5))) = -3;
    ipHat(find(r>= 2/sqrt(5))) = 3;
    ipHat(find(r>=-2/sqrt(5) & r<0)) = -1;
    ipHat(find(r>=0 & r<2/sqrt(5))) = 1;
    nErr(i) = size(find([ip- ipHat]),2); % couting the number of errors
end
simBer = nErr/N;

% plot
figure
hold on
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
axis([0 10 10^-4 1])
grid on
legend('zf-5tap','zf-7tap','zf-9tap','mmse-5tap','mmse-7tap','mmse-9tap','AWGN');
xlabel('Eb/No, dB');
ylabel('Bit Error Rate');
title('Bit error probability curve for 4-PAM in ISI with ZF/MMSE equalizer');