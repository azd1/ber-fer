% ECE 6606 Project 2
% MESONERO Philippe 
% ZDUN Arnaud

% NB:
% PE   = power efficiency 
% PEdB = power efficiency in Decibel
% rho  = spectral efficiency
clear;
close all;
clc;
%% 1  Benchmarks: Shannon limit andQPSK limit

max = 12;
step = 1e-3;
dim = max/step + 1;

% Equation Number 1
rho = 0:step:max;

% Equation number 2 

PEdB2 = zeros(1,dim);
N0 = 2*10.^(-2:0.001:1);
sigma2 = N0./2;


% Trade-Off Equation 1
PE1 = (2.^rho - ones(1,dim)) ./ rho;
PE1dB = 10 * log10(PE1);

% Calculation for Eq 2 
p = @(y) (exp(-(y.^2 + ones(1,length(y))) ./ (2 * sigma2)) .* cosh(y./sigma2) ./ sqrt(2 * pi .* sigma2));
fsigma2 = -2 * integral(@(y) (p(y) .* log2(p(y))),-8,8,'ArrayValued',true) - log2(2 * pi * sigma2 * exp(1));

PE2 = ones(1,length(sigma2)) ./ (sigma2 .* fsigma2);
PE2dB = 10 * log10(PE2);

index = rho>=0 & rho<=10 & PE1dB>=-2 & PE1dB<=12;  
index2 = fsigma2>=0.5 & fsigma2<=3 & PE2dB>=-2 & PE2dB<=12;  

figure
plot (PE1dB(index), rho(index), PE2dB(index2), fsigma2(index2), 10.4, 7/8, '+', 9.5, 15/22, '+',9.1, 31/52, '+', 11.3, 1, 'd', 9.6, 2/3, 'v',10.1, 2/11, 'v', 'LineWidth', 1.5);
grid on
xlabel('Power Efficiency $E_b/N_0$ (dB)','interpreter','latex')
ylabel('Spectral Efficiency $\rho$ (b/s/Hz)','interpreter','latex')
l = legend('Fundamental Trade off without any constraint','Fundamental Trade off for QPSK','Hamming (7,4)','Hamming (15,11)','Hamming (31,26)','Golay','3-Repetition','11-Repetition');
l.FontName = 'Times New Roman';

%% Theorical BER Uncoded QPSK

% Theorical BER of an uncoded QPSK
q = @(x) integral(@(t)(exp(-t.^2 / 2) / sqrt(2 * pi)), x, 10);

%index3 = PE2dB>=-2 & PE2dB<=15;

% Theoritical BER
BER = arrayfun(q, sqrt(0.5 ./ sigma2));
BERdB = 10*log10(BER);
PowEff0 = PowEff(0.5);
figure
loglog(PowEff0,BER)
grid on
title ('Theoritical Bit Error Rate for uncoded QPSK')
xlabel('Power Efficiency Eb/No(dB)')
ylabel('BER')

%%
close all;
sigma2 = 10.^(-2:0.01:1);
PowEff = @(EB) EB ./ (2*sigma2);

%% uncoded
n = 8;
k = 8;
Es = 1; % energy per QPSK symbol
Eb = Es/2*(n/k);
poweffunc = PowEff(Eb);
[ferRunc,berRunc] = computeBFER(n, k, Es, Eb, sigma2, 40000, 'none', 'berferUnc.mat', true);

%% repetition codes
n = 3;
k = 1;
Es = 1; % energy per QPSK symbol
Eb = n*Es/2; % energy per coded bit
poweff1 = PowEff(Eb);
[ferR3,berR3] = computeBFER(n, k, Es, Eb, sigma2, 10000, 'repetition', 'berferRepet.mat', true);

n = 11;
k = 1;
Es = 1; % energy per QPSK symbol
Eb = n*Es/2; % energy per coded bit
poweff2 = PowEff(Eb);
[ferR11,berR11] = computeBFER(n, k, Es, Eb, sigma2, 10000, 'repetition', 'berferRepet11.mat', true);
%% Plot part 1 
figure
semilogy(10*log10(poweff1),berR3, 10*log10(poweff2),berR11);
grid on
%title('Empirical Bit Error Rate of the repetion code [n,1,n],  respect to the power efficiency','interpreter','latex')
xlabel ('$E_b/N_0$ (dB)','interpreter','latex')
ylabel('Bit Error Rate','interpreter','latex')
l = legend('n = 3','n = 11');
l.FontName = 'Times New Roman';


figure
semilogy(10*log10(poweff1),ferR3, 10*log10(poweff2),ferR11);
grid on
%title('Empirical Frame Error Rate of the repetion code [n,1,n],  respect to the power efficiency','interpreter','latex')
xlabel ('$E_b/N_0$ (dB)','interpreter','latex')
ylabel('Frame Error Rate','interpreter','latex')
l = legend('n = 3','n = 11');
l.FontName = 'Times New Roman';

%% Hamming
%compute Hamming (n,k)
n = 7;
k = 4;
Es = 1; % energy per QPSK symbol
Eb = n*Es/(2*k); % energy per coded bit
poweffH74 = PowEff(Eb);
[ferH74,berH74] = computeBFER(n, k, Es, Eb, sigma2, 1, 'Hamming', 'Hamming74.mat', true);

n = 15;
k = 11;
Es = 1; % energy per QPSK symbol
Eb = n*Es/(2*k); % energy per coded bit
poweffH1511 = PowEff(Eb);
[ferH1511,berH1511] = computeBFER(n, k, Es, Eb, sigma2, 1, 'Hamming','Hamming1511.mat', true);

n = 31;
k = 26;
Es = 1; % energy per QPSK symbol
Eb = n*Es/(2*k); % energy per coded bit
poweffH3126 = PowEff(Eb);
[ferH3126,berH3126] = computeBFER(n, k, Es, Eb, sigma2, 1, 'Hamming','Hamming3126.mat', true);

%%
n = 24;
k = 12;
Es = 1; % energy per QPSK symbol
Eb = n*Es/(2*k); % energy per coded bit
poweffGolay = PowEff(Eb);
[ferGolay,berGolay] = computeBFER(n, k, Es, Eb, sigma2, 1000, 'Golay','Golay3.mat', true);

%%
q = @(x) integral(@(t)(exp(-t.^2 / 2) / sqrt(2 * pi)), x, 10);
% Theoretical BER
BER = arrayfun(q, sqrt(0.5 ./ sigma2));
poweff = PowEff(0.5);
le = length(poweff); step = ceil(le/100);
%semilogy(10*log10(poweffH74),berH74, 10*log10(poweffH1511),berH1511, 10*log10(poweffH3126),berH3126, 10*log10(poweff(1:step:end)),BER(1:step:end),'o--','LineWidth',1.5);
semilogy(10*log10(poweffGolay), berGolay, 10*log10(poweff(1:step:end)),BER(1:step:end),'+--','LineWidth',1.5);
%semilogy(10*log10(poweffunc),berRunc, 10*log10(poweff(1:step:end)),BER(1:step:end),'o--','LineWidth',1.5);
grid on
%title('Empirical Frame Error Rate of the repetion code [n,1,n],  respect to the power efficiency','interpreter','latex')
xlabel ('$E_b/N_0$ (dB)','interpreter','latex')
ylabel('Frame Error Rate','interpreter','latex')
l = legend('Golay code', 'Uncoded');
%l = legend('H(7,4)', 'H(15,11)','H(31,26)', 'Uncoded QPSK');
l.FontName = 'Times New Roman';
axis([-2,12,1E-5,1]);

%%
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'interval','-dpdf','-r0')

