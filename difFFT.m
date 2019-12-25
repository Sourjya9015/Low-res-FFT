% Do a limited precision 2048/1024 point DIF FFT
clc; clear;
nfft = 2048;
nsym = 1;
inPrec = 4; % 4 bits of precision on the input
inFrac = inPrec - 2; % 2 bits for fraction, 1 bit integer and 1 bit sign

finalPrec = 8;
finalFrac = 4;

% F1 = fimath ('RoundingMethod','Nearest', 'OverflowAction', 'Saturate', ...
%                  'ProductMode', 'FullPrecision' , 'SumMode', 'FullPrecision');
             
F1 = fimath ('RoundingMethod','Nearest', 'OverflowAction', 'Saturate', ...
                 'ProductMode', 'SpecifyPrecision' , ...
                 'ProductFractionLength', 10,...
                 'ProductWordLength', 16, ...
                 'SumMode', 'SpecifyPrecision',...
                 'SumFractionLength', 10,...
                 'SumWordLength', 16);

% Assume input data to be between 0 and 1
% reduces the error due to the initial inPrec bit 
% truncation of the data.

% random and normalized signal
% xti = randn(nfft,nsym) + 1i*randn(nfft,nsym);
% xmax = max(max(xti));
% xti = xti./xmax;

% Sine
% t = (1:(nfft*nsym))';
% w0 = 2*pi*0.3;  w1 = 2*pi*0.7;
% xti = 0.25*sin(w0*t) + 0.5*sin(w1*t) + 1i*0.25*cos(w0*t) + 1i*0.5*cos(w1*t);
% xti = xti*8;

% QPSK
xfreqDom = 16*exp(1i*pi/2*(randi(4,nfft,nsym)+0.5));
xti = ifft(xfreqDom, nfft);



% Permute the data
indx = 0:(nfft-1);
binIndx = dec2bin(indx);
binIndx = fliplr(binIndx);
indx = bin2dec(binIndx) + 1;
%xti = xti(indx,:);

xtf = fi(xti,1,inPrec, inFrac, F1);

Xf = fft(double(xtf), nfft);
xt = xtf(indx,:);

% Stage 1 and 2 have only additions and trivial phase multiplications
% i.e. with 1+0i, -1+0i, 0+1i, 0-1i.

% stage 1: No complex multipliers required!
xt1 = zeros(nfft,1) + 1i*zeros(nfft,1);
xt2 = zeros(nfft,1) + 1i*zeros(nfft,1);

xt1 = fi(xt1, 1, inPrec+1, inFrac,F1);
for i=1:2:nfft
    temp = xt(i+1);
    xt1(i+1) = xt(i) - temp;
    xt1(i) = xt(i) + temp;
end
%disp('Stage 1');
%disp(xt1.WordLength); disp(xt1.FractionLength);
% Stage 2
w2 = exp(-1i*2*pi*(0:(2^2/2-1))/(2^2));
w2 = fi(w2,1,2,0,F1); % two bits are enough in this case

xt2 = fi(xt2, 1, xt1.WordLength+1, xt1.FractionLength,F1);
%for mcnt=1:2
    for j=1:4:nfft
        indx = j+2;
        temp2 = xt1(indx);
        xt2(indx) = xt1(j) - temp2;
        xt2(j) = xt1(j) + temp2;
    end
    
    for j=2:4:nfft
        indx = j+2;
        temp2 = imag(xt1(indx)) - 1i*real(xt1(indx));
        xt2(indx) = xt1(j) - temp2;
        xt2(j) = xt1(j) + temp2;
    end
%end

%disp('Stage 2');
%disp(xt2.WordLength); disp(xt2.FractionLength);

M = 4;
             
wordLen = 0;
fracLen = 0;

while M < nfft
    wM = exp(-1i*2*pi*(0:(M-1))/(2*M));
    %wM = fi(wM,1,8,6,F1);
    if M < nfft/2
        wM = fi(wM,1,10,8,F1);
    else
        wM = fi(wM,1,8,6,F1);
    end
    
    Istep = M*2;
    xt3 = zeros(nfft,1) + 1i*zeros(nfft,1);
    
    if ( xt2.WordLength < finalPrec-8)
        wordLen = xt2.WordLength+4;
        fracLen = xt2.FractionLength+2;
    else
        wordLen = finalPrec;
        fracLen = finalFrac;
    end

    xt3 = fi(xt3, 1, wordLen, fracLen,F1);
    for mcnt = 1:M
        for icnt = mcnt:Istep:nfft 
            j = icnt + M;
            Temp = wM(mcnt)* xt2(j);
            xt3(j) = xt2(icnt) - Temp;
            xt3(icnt) = xt2(icnt) + Temp;
        end
    end
    xt2 = xt3;
    M = Istep;
    %disp('Stage 2');
    %disp(xt2.WordLength); disp(xt2.FractionLength);
end

xt2 = fi(xt2,1,4);
figure(9);
plot(real(xt2), imag(xt2),'x'); hold on;
plot(real(Xf), imag(Xf),'or');

disp(mean(abs(Xf-double(xt2))));

xfLowPrec = double(xt2);

% 
% xfLowPrec = fi(xfLowPrec, 1,6);
% xfLowPrec = double(xfLowPrec);
% % Plot the FFT spectrums
% powXf = abs(Xf);
% powXf0 = abs(xfLowPrec);
% 
% nbins = length(powXf0);
% semilogy(0:1/(nbins/2-1):1, powXf(1:nbins/2),'-r'); hold on;
% semilogy(0:1/(nbins/2-1):1, powXf0(1:nbins/2),'-.');
% legend('\infty Precision','Finite precision');
% set(gca,'Fontsize',20);
% xlabel('\omega/pi'); ylabel('|P(\omega)|');
% grid on;

% for QPSK
phaseRx = angle(xfLowPrec);
phaseTx = angle(xfreqDom);
phaseHighPrec = angle(Xf);

phaseDiff = abs(phaseTx - phaseRx);
phaseDiff4 = abs(phaseTx - phaseHighPrec);

errTot = sum(phaseDiff > pi/4);
errTot4 = sum(phaseDiff4 > pi/4);
%err = errTot/nfft/nsym;
fprintf('Number of symbols in error = %d\n',errTot);
fprintf('Symbols in error due to 4 bit truncation = %d\n',errTot4);





