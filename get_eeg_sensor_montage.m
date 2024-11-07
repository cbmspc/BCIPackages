% Master channels list
%
%

function [Montage, MontageDesc, NexusNumS] = get_eeg_sensor_montage ()
Montage.nexus1 = string_to_cell('Cz, Fp1, Fpz, Fp2, F7, F3, Fz, F4, F8, FC5, FC1, FC2, FC6, M1, T7, C3, C4, T8, M2, CP5, CP1, CP2, CP6, P7, P3, Pz, P4, P8, POz, O1, Oz, O2, D1', ', ').';
MontageDesc.nexus1 = 'Nexus Amp #1 default';
NexusNumS.nexus1 = 1;

Montage.nexus2 = string_to_cell('Cz, AF7, AF3, AFz, AF4, AF8, F5, F1, F2, F6, FC3, FCz, FC4, C5, C1, C2, C6, CP3, CPz, CP4, P5, P1, P2, P6, PO3, PO4, FT7, FT8, TP7, TP8, PO7, PO8, D2', ', ').';
MontageDesc.nexus2 = 'Nexus Amp #2 default';
NexusNumS.nexus2 = 2;

Montage.numericnexus1 = string_to_cell([num2str(1:32,'A%i ') ' D1'],' ').';
MontageDesc.numericnexus1 = 'Nexus Amp #1 numeric';
NexusNumS.numericnexus1 = 6;

Montage.numericnexus2 = string_to_cell([num2str(1:32,'B%i ') ' D2'],' ').';
MontageDesc.numericnexus2 = 'Nexus Amp #2 numeric';
NexusNumS.numericnexus2 = 10;

Montage.synfi = [Montage.nexus1; Montage.nexus2; {'DIAG'}];
MontageDesc.synfi = 'SynFi default';
NexusNumS.synfi = 3;

Montage.numericsynfi = [Montage.numericnexus1; Montage.numericnexus2; {'DIAG'}];
MontageDesc.numericsynfi = 'SynFi numeric';
NexusNumS.numericsynfi = 7;

Montage.nexus1 = [Montage.nexus1; {'DIAG'}];
Montage.numericnexus1 = [Montage.numericnexus1; {'DIAG'}];

Montage.nexus2 = [Montage.nexus2; {'DIAG'}];
Montage.numericnexus2 = [Montage.numericnexus2; {'DIAG'}];

Montage.nexusbbs = string_to_cell('FP1, FPZ, FP2, AF7, AF3, AFZ, AF4, AF8, F7, F5, F3, F1, FZ, F2, F4, F6, F8, FT7, FC5, FC3, FC1, FCZ, FC2, FC4, FC6, FT8, T7, C5, C3, C1, CZ, C2, C4, C6, T8, TP7, CP5, CP3, CP1, CPZ, CP2, CP4, CP6, TP8, P7, P5, P3, P1, PZ, P2, P4, P6, P8, PO7, PO3, POZ, PO4, PO8, O1, OZ, O2, M1, M2', ', ').';
MontageDesc.nexusbbs = 'BioTrace default';
NexusNumS.nexusbbs = 11;

Montage.zoranoc1 = string_to_cell('AFz, F3, F1, Fz, F2, F4, FC5, FC3, FC1, FCz, FC2, FC4, FC6, C5, C3, C1, Cz, C2, C4, C6, CP5, CP3, CP1, CPz, CP2, CP4, CP6, P3, P1, Pz, P2, P4, D1, DIAG', ', ').';
MontageDesc.zoranoc1 = 'Zoran OC1';
NexusNumS.zoranoc1 = 4;

Montage.zoranoc2 = string_to_cell('C1, C2, C3, C4, C5, C6, CP1, CP2, CP3, CP4, CP5, CP6, CPz, Cz, F1, F2, F3, F4, FC1, FC2, FC3, FC4, FC5, FC6, FCz, Fz, P1, P2, P3, P4, POz, Pz, D1, DIAG', ', ').';
MontageDesc.zoranoc2 = 'Zoran OC2';
NexusNumS.zoranoc2 = 5;

Montage.acticap1 = string_to_cell('FP1 FZ F3 F7 FT9 FC5 FC1 C3 T7 M1 CP5 CP1 PZ P3 P7 O1 OZ O2 P4 P8 M2 CP6 CP2 CZ C4 T8 FT10 FC6 FC2 F4 F8 FP2 D1 DIAG', ' ').';
MontageDesc.acticap1 = 'ActiCap1';
NexusNumS.acticap1 = 9;

Montage.zoranoc3 = string_to_cell('F5, F1, Fz, F2, F6, FC3, FCz, FC4, C5, C3, C1, Cz, C2, C4, C6, CP3, CPz, CP4, P5, P1, Pz, P2, P6, POz, NC1-25, NC1-26, NC1-27, NC1-28, NC1-29, NC1-30, NC1-31, NC1-32, D1, DIAG', ', ').';
MontageDesc.zoranoc3 = 'Zoran OC3';
NexusNumS.zoranoc3 = 12;

Montage.brainlesion = string_to_cell('Cz, Fz, Pz, Oz, P3, P4, O1, O2, M1, M2, C3, C4, NC1-13, NC1-14, NC1-15, NC1-16, NC1-17, NC1-18, NC1-19, NC1-20, NC1-21, NC1-22, NC1-23, NC1-24, NC1-25, NC1-26, NC1-27, NC1-28, NC1-29, NC1-30, NC1-31, NC1-32, D1, DIAG', ', ').';
MontageDesc.brainlesion = 'brainlesion';
NexusNumS.brainlesion = 13;

ChansGreen = {'Fp1','Fp2','F7','F3','Fz','F4','F8','FC5','FC1','FC2','FC6','T7','C3','Cz','C4','T8','TP9','CP5','CP1','CP2','CP6','TP10','P7','P3','Pz','P4','P8','PO9','O1','Oz','O2','PO10'};
ChansRed = {'Fpz','F9','AFF5h','AFF1h','AFF2h','AFF6h','F10','FTT9h','FTT7h','FCC5h','FCC3h','FCC1h','FCC2h','FCC4h','FCC6h','FTT8h','FTT10h','TPP9h','TPP7h','CPP5h','CPP3h','CPP1h','CPP2h','CPP4h','CPP6h','TPP8h','TPP10h','POO9h','POO1','POO2','POO10h','Iz'};
ChansYellow = {'AF7','AF3','AF4','AF8','F5','F1','F2','F6','FT9','FT7','FC3','FC4','FT8','FT10','C5','C1','C2','C6','TP7','CP3','CPz','CP4','TP8','P5','P1','P2','P6','PO7','PO3','POz','PO4','PO8'};
ChansWhite = {'AFp1','AFp2','FFT9h','FFT7h','FFC5h','FFC3h','FFC1h','FFC2h','FFC4h','FFC6h','FFT8h','FFT10h','TTP7h','CCP5h','CCP3h','CCP1h','CCP2h','CCP4h','CCP6h','TTP8h','P9','PPO9h','PPO5h','PPO1h','PPO2h','PPO6h','PPO10h','P10','I1','OI1h','OI2h','I2'};
Montage.acticap128_1 = [ChansGreen,ChansYellow,ChansRed,ChansWhite].';
MontageDesc.acticap128_1 = 'ActiCap128_1';
NexusNumS.acticap128_1 = 14;

Montage.waveguard = string_to_cell('Cz, Fp1, Fpz, Fp2, F7, F3, Fz, F4, F8, FC5, FC1, FC2, FC6, M1, T7, C3, C4, T8, M2, CP5, CP1, CP2, CP6, P7, P3, Pz, P4, P8, POz, O1, Oz, O2, TRIG_1, NC, AF7, AF3, NC, AF4, AF8, F5, F1, F2, F6, FC3, NC, FC4, C5, C1, C2, C6, CP3, CPz, CP4, P5, P1, P2, P6, PO3, PO4, FT7, FT8, TP7, TP8, PO7, PO8, TRIG_2, DIAG', ', ').';
MontageDesc.waveguard = 'waveguard';
NexusNumS.waveguard = 15;

NexusNumS.userdefined = 8;

Montage.MCN = {'Fp1', 'Fpz', 'Fp2', 'AF7', 'AF3', 'AFz', 'AF4', 'AF8', 'F9', 'F7', 'F5', 'F3', 'F1', 'Fz', 'F2', 'F4', 'F6', 'F8', 'F10', 'FT9', 'FT7', 'FC5', 'FC3', 'FC1', 'FCz', 'FC2', 'FC4', 'FC6', 'FT8', 'FT10', 'T7', 'C5', 'C3', 'C1', 'Cz', 'C2', 'C4', 'C6', 'T8', 'TP9', 'TP7', 'CP5', 'CP3', 'CP1', 'CPz', 'CP2', 'CP4', 'CP6', 'TP8', 'TP10', 'P9', 'P7', 'P5', 'P3', 'P1', 'Pz', 'P2', 'P4', 'P6', 'P8', 'P10', 'PO7', 'PO3', 'POz', 'PO4', 'PO8', 'O1', 'Oz', 'O2'};
MontageDesc.MCN = 'Modified Combinatorial Nomenclature';
NexusNumS.MCN = 16;

