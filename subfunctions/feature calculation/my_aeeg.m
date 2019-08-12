function [t, aeeg, aeeg_nl, aref] = my_aeeg(eeg, fs, flag)
% A function for estimating the amplitude integrated EEG
% Two estimates are provided and selected using the flag variable.
%
% INPUTS: eeg - a vector containing one-channel of EEG 
%                  fs - a scalar containing the sampling frequency of the EEG
%                  flag - a scalar that selectes the type of aEEG estimate
%                  so that 0 provides a crude approximation to a CFM
%                  machine, and 1, provides an approximation to the Olympic CFM monitor
%
%  OUTPUTS: t - a vector containt the time
%                      aeeg - a vector containing the aEEG value if flag is
%                      0 then this is one row vector containing the
%                      CFM output and if flage is 1 then a two-row vector
%                      is output containt the upper and lower margins of
%                      the aEEG estimate.
%                      aeeg_nl - the nonlinearly scaled version of the aeeg for display
%                      aref - reference values to correctly label plots.
%
%
% NOTE that the aEEG variable from method flag=0 will be most similar to
% the dum1 variable in method flag~=0. dum1 is not passed from this
% method, however, as quantiles are passed instead. Feel free to modify
% code accordingly.
%
%
% Nathan Stevenson
% University of Helsinki, Finland
% 15/03/02017


if flag==0
    
    % OLD SCHOOL CFM VERSION

load aeeg_filter % EEG filter from CFM machine pass band implemented with fs=256
% use freqz(b,a,4096, 256) to see frequency response sampling frequency=256
eeg = resample(eeg, 256, fs); % resample EEG to fit filter coefficients
ef = filter(b,a,eeg);
epl = 1; olap = 0.5; % variables relating to average time and temporal resoultion of aEEG estimate
block_no = floor(length(ef)/(olap*256)-(epl/olap));
aeeg = zeros(1,block_no);
for ii = 1:block_no;
    r1 = (ii-1)*olap*256+1; r2 = r1+epl*256-1;
    aeeg(ii) = 2*mean(abs(ef(r1:r2)));  % *2 to generate a peak to peak estimate
end

% nonlinear mapping function in old CFMs it was not logarithmic rather
% piecewise linear
aeeg_nl = nlin_map(tp1, tp2, s1, s2, s3, aeeg);
t = 0:olap:length(aeeg)*olap-olap; % time vector

% nonlinear mapping to help scale axis correctly
ef1 = [0 5 10 20 50 100];
aref(1,:) = ef1;
aref(2,:) = nlin_map(tp1, tp2, s1, s2, s3, ef1);

else

    % DESIGNED TO APPROXIMATE OLYMPIC CFM

    % aEEG fiter - similar to above but FIR and will work for an arbitray
    % sampling frequency
Ny = fs/2;
Fo =  [0 .05 2 15 21 Ny]/Ny;
Mo =  [0 0 0.3 1.2 0 0];
W =   [30 .3 30];
B = firpm(Ny,Fo,Mo,W);  % Parks-McCellan optimation FIR filter design
yF = filter(B,1,eeg);
yF1 = 2*abs(yF); % *2 for peak to peak estimate

epl = 3; 
% smooth estimate using 3s rectanguler window -  this is basically the aEEG
dum1 = conv(yF1, ones(1,epl*fs));
dum1 = dum1(epl*fs/2:end-epl*fs/2)./(epl*fs);

% this stage approximates the process seen in some machines like the NicOne
% which plot a lower and upper quantile rather than the aEEG.
epl = 30; olap = 1;
block_no2 = floor(length(dum1)/(olap*fs)-(epl/olap));
aeeg = zeros(2,block_no2);
for ii = 1:block_no2
    r1 = (ii-1)*olap*fs+1; r2 = r1+epl*fs-1;
    aeeg(1,ii) = quantile(dum1(r1:r2), 0.1); 
    aeeg(2,ii) = quantile(dum1(r1:r2), 0.9); 
end
t = 0:olap:length(aeeg)*olap-olap;
aeeg_nl = aeeg;
% Nonlinear mapping, this is the more standard semi-log function that is
% linear<10 then log>10
aeeg_nl(1,aeeg(1,:)>10) = 10*log10(aeeg(1,aeeg(1,:)>10));
aeeg_nl(2,aeeg(2,:)>10) = 10*log10(aeeg(2,aeeg(2,:)>10));

ef1 = [0 5 10 20 50 100];
aref(1,:) = ef1; aref(2,:) = ef1;
aref(2, ef1>10) = 10*log10(ef1(ef1>10));


end

