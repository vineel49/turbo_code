% Turbo coded QPSK OVER AWGN - overall rate is 1/2

close all
clear all
clc
SNR_dB = 3; % SNR PER BIT in dB
FRAME_SIZE = 1024;
NUM_BIT = 0.5*FRAME_SIZE; % NUMBER OF DATA BITS - OVERALL RATE IS 1/2
NUM_FRAMES = 10^1; % NUMBER OF FRAMES SIMULATED

% SNR PARAMETERS - OVERALL RATE IS 1/2
SNR = 10^(0.1*SNR_dB); % SNR IN LINEAR SCALE
NOISE_VAR_1D = 2*2/(2*SNR); % 1D AWGN NOISE VARIANCE 
NOISE_STD_DEV = sqrt(NOISE_VAR_1D); % NOISE STANDARD DEVIATION
%--------------------------------------------------------------------------
% GENERATOR polynomial of the component encoders
GEN_POLY = ldiv2([1 0 1],[1 1 1],NUM_BIT); % using long division method

% Interleaver and deinterleaver mapping of the turbo code 
INTR_MAP = randperm(NUM_BIT);
DEINTR_MAP = deintrlv((1:NUM_BIT),INTR_MAP);

tic()
C_BER = 0; % bit errors in each frame
for FRAME_CNT = 1:NUM_FRAMES
%----            TRANSMITTER      -----------------------------------------
% SOURCE
A = randi([0 1],1,NUM_BIT);

% Turbo encoder
% component encoder 1
B1 = zeros(1,2*NUM_BIT); % encoder 1 output initialization
B1(1:2:end) = A; % systematic bit
temp1 = mod(conv(GEN_POLY,A),2); 
B1(2:2:end) = temp1(1:NUM_BIT); % parity bit
% component encoder 2
B2 = zeros(1,2*NUM_BIT); % encoder 2 output initialization
B2(1:2:end) = A(INTR_MAP); % systematic bit
temp2 = mod(conv(GEN_POLY,B2(1:2:end)),2); 
B2(2:2:end) = temp2(1:NUM_BIT); % parity bit

% QPSK MAPPING
% QPSK mapping (according to the set partitioning principles)
MOD_SIG1 = 1-2*B1(1:2:end) + 1i*(1-2*B1(2:2:end));
MOD_SIG2 = 1-2*B2(1:2:end) + 1i*(1-2*B2(2:2:end));
MOD_SIG = [MOD_SIG1 MOD_SIG2];

%---------------     CHANNEL      -----------------------------------------
% AWGN
AWGN = normrnd(0,NOISE_STD_DEV,1,FRAME_SIZE)+1i*normrnd(0,NOISE_STD_DEV,1,FRAME_SIZE);

% CHANNEL OUTPUT
CHAN_OP = MOD_SIG + AWGN;
%----------------      RECEIVER  ------------------------------------------
%--------------------ITERATIVE TURBO DECODING -----------------------------
% Branch metrices for the SOVA
QPSK_SYM = zeros(4,2*NUM_BIT);
QPSK_SYM(1,:) = (1+1i)*ones(1,2*NUM_BIT);
QPSK_SYM(2,:) = (1-1i)*ones(1,2*NUM_BIT);
QPSK_SYM(3,:) = (-1+1i)*ones(1,2*NUM_BIT);
QPSK_SYM(4,:) = (-1-1i)*ones(1,2*NUM_BIT);

Dist = zeros(4,2*NUM_BIT);
 Dist(1,:)=abs(CHAN_OP-QPSK_SYM(1,:)).^2;
 Dist(2,:)=abs(CHAN_OP-QPSK_SYM(2,:)).^2;
 Dist(3,:)=abs(CHAN_OP-QPSK_SYM(3,:)).^2;
 Dist(4,:)=abs(CHAN_OP-QPSK_SYM(4,:)).^2;
  
BRANCH_METRIC1 = Dist(:,1:NUM_BIT); % branch metrices for component decoder 1
BRANCH_METRIC2 = Dist(:,NUM_BIT+1:end); % branch metrices for component decoder 2
 
 % a priori probabilities (LLR) - initialization
APR_LLR = zeros(1,NUM_BIT); % for first iteration

% iterative decoding
SOFT_OUTPUT = SOVA(APR_LLR,NUM_BIT,BRANCH_METRIC1);
SOFT_OUTPUT = SOVA(SOFT_OUTPUT(INTR_MAP),NUM_BIT,BRANCH_METRIC2); %1

SOFT_OUTPUT = SOVA(SOFT_OUTPUT(DEINTR_MAP),NUM_BIT,BRANCH_METRIC1);
SOFT_OUTPUT = SOVA(SOFT_OUTPUT(INTR_MAP),NUM_BIT,BRANCH_METRIC2); %2

SOFT_OUTPUT = SOVA(SOFT_OUTPUT(DEINTR_MAP),NUM_BIT,BRANCH_METRIC1);
SOFT_OUTPUT = SOVA(SOFT_OUTPUT(INTR_MAP),NUM_BIT,BRANCH_METRIC2); %3

SOFT_OUTPUT = SOVA(SOFT_OUTPUT(DEINTR_MAP),NUM_BIT,BRANCH_METRIC1);
SOFT_OUTPUT = SOVA(SOFT_OUTPUT(INTR_MAP),NUM_BIT,BRANCH_METRIC2); %4

SOFT_OUTPUT = SOVA(SOFT_OUTPUT(DEINTR_MAP),NUM_BIT,BRANCH_METRIC1);
SOFT_OUTPUT = SOVA(SOFT_OUTPUT(INTR_MAP),NUM_BIT,BRANCH_METRIC2); %5

SOFT_OUTPUT = SOVA(SOFT_OUTPUT(DEINTR_MAP),NUM_BIT,BRANCH_METRIC1);
SOFT_OUTPUT = SOVA(SOFT_OUTPUT(INTR_MAP),NUM_BIT,BRANCH_METRIC2); %6

SOFT_OUTPUT = SOVA(SOFT_OUTPUT(DEINTR_MAP),NUM_BIT,BRANCH_METRIC1);
SOFT_OUTPUT = SOVA(SOFT_OUTPUT(INTR_MAP),NUM_BIT,BRANCH_METRIC2); %7

SOFT_OUTPUT = SOVA(SOFT_OUTPUT(DEINTR_MAP),NUM_BIT,BRANCH_METRIC1);
SOFT_OUTPUT = SOVA(SOFT_OUTPUT(INTR_MAP),NUM_BIT,BRANCH_METRIC2); %8

SOFT_OUTPUT = SOVA(SOFT_OUTPUT(DEINTR_MAP),NUM_BIT,BRANCH_METRIC1);
SOFT_OUTPUT = SOVA(SOFT_OUTPUT(INTR_MAP),NUM_BIT,BRANCH_METRIC2); %9

SOFT_OUTPUT = SOVA(SOFT_OUTPUT(DEINTR_MAP),NUM_BIT,BRANCH_METRIC1);
SOFT_OUTPUT = SOVA_END(SOFT_OUTPUT(INTR_MAP),NUM_BIT,BRANCH_METRIC2); %10

% hard decision is taken on the a posteriori LLR
SOFT_OUTPUT  = SOFT_OUTPUT(DEINTR_MAP);
DEC_A = SOFT_OUTPUT<0;
% CALCULATING BIT ERRORS IN EACH FRAME
C_BER = C_BER + nnz(A-DEC_A);
end
toc()
% bit error rate
BER = C_BER/(NUM_BIT*NUM_FRAMES)