clc;clear;
nfft = 2048;
w = {};
pos = 1;
M = 1;
while M < nfft
    w{pos} = exp(-1i*2*pi*(0:(M-1))/(2*M));
    M = M*2;
    pos = pos+1;
end

% w is a 11 element cell array. Each cell is a variable length vector.