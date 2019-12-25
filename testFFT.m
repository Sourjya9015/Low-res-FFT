
clc; clear;


nprecArr = 10;% 10 12 16];
nfft = [16 64 256 512 1024 2048];

err = zeros(length(nprecArr),length(nfft));

for j=1:length(nprecArr)
    nprec = nprecArr(j);
    for icnt=1:length(nfft)
        At = randn(nfft(icnt),5) + 1i*randn(nfft(icnt),5);
        
        Atfin = fi(At, 1, 4, 2);
        At = double(Atfin);

        trueFFT = fft (At, nfft(icnt));

        lowPrecFFt = computeLowPrecFFT( At, nfft(icnt), nprec);

        err(j,icnt) = mean(mean((abs(trueFFT - double(lowPrecFFt))).^2));

    end
end

figure(1);
semilogy(nfft, err, '-v', 'Linewidth', 2);
xlabel('N_{FFT}'); ylabel('MSE');
%legend('n=8','n=10','n=12','n=16');
grid on;


