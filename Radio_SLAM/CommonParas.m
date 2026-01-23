% System Parameters---------------------
UENum = 16;
LightSpeed = 3e8;
SNR1_dB = 20;
SNR2_dB = 10;
SubcarrierNum = 2048;
SubcarrierSpacing = 30e3;

% Method Parameters---------------------
% Extract path---
delta_OMP = 1e-6;
MaxIterNum_OMP = 15;
MaxDelaySpread = 200e-9;
OverSamplingFactor = 1;
% SLAM
MinANum = 5;
SolveWallTols = [1,1,1];
PowerCheckThres = 10;%dB
RandNumPerSeed = 300;
