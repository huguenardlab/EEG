% ----- [fd,param_fund,fit_fund,err] = harmonic_filter(d,Fs,varargin) -----
%
%
%   Harmonic Filter used to cancel out noise and its harmonics.  Typical
%   filters normally implement bandpasses and bandstops to filter out noise
%   from a signal, however this can lead to artifacts and fluctuations at
%   deflection points in the data.  Instead, harmonic filters fit
%   sinusoiuds to the fundamental frequency and harmonics of the 
%   confounding noise and subtract them from the original signal.
%
%   This function makes use of "sinefit.m", which fits sinusoids without
%   specific parameters (via least squared error).  
%   I wonder if it might do better generally, if the signal was first high pass filtered
% at 1/2 or so of the fundamental?  JRH 11/21/18 
% then even if there were drift in the baseline, the sine fit would still settle on the noise, 
% especially if it were phase and amplitude invariant.  The latter would need to be true generally in order 
% for a subtraction of a single sine wave to work.
%
%                   >>> INPUTS >>>
%   REQUIRED:
%       d = data (vector or array)
%       Fs = sampling rate
%   OPTIONAL:
%       fund = fundamental frequency of noise (default = 60 Hz)
%       harm = # of harmonics to fit and subtract (default = 1)
%       plotting = 0 or 1 (default 0)
%           - plots original and filtered signal if 1
%
%                   <<< OUTPUTS <<<
%       fd = de-noised data
%       params = fitting parameters
%       fit = actual sinusoid fit
%       err = LSE value
%
%   Examples:
%   (1) [fd] = harmonic_filter(d,Fs,[],2);
%           - filters out 60hz and 2 harmonics
%
%   (2) [fd] = harmonic_filter(d,Fs,45,3,1);
%           - filters out 45hz and 3 harmonics, and plots the results
%       
%
% By JMS, 2/17/2015
% -----------------------------------------------------------------------

function [fd] = harmonic_filter2(d,Fs,varargin)

    % extract optional inputs
    if nargin>4;plotting=varargin{3};
    else plotting=0;end
    if nargin>3;harmonics=varargin{2};
    else harmonics=0;end
    if nargin>2;fundamental=varargin{1};
    else fundamental=60;end

    % error checking
    if isempty(fundamental);fundamental=60;end
    if isempty(harmonics);harmonics=1;end
    if isempty(plotting);plotting=0;end
    if fundamental<0;fundamental=60;end
    if harmonics<0;harmonics=0;end;
    if plotting<0 || plotting>1;plotting=0;end

    % preprocessing
    si = 1/Fs;
    stop = length(d)/Fs; % amount of data to be filtered
    windowlength=length(d);
    
    maxfilterinputlength=1e3;
       if length(d)>maxfilterinputlength
        windowlength=maxfilterinputlength;    
       end
    windows=floor(length(d)/windowlength);
    time=(0:windowlength-1)/Fs;
    fdout=zeros(windows*windowlength,1);
    fund = fundamental-1:0.5:fundamental+1; % +/- 2 from fundamental frequency
    
    %fprintf('\nRemoving hum from %d windows, window #%5d',windows,0);
f = waitbar(0,'Removing line hum');
winupdate=0;
    for window=1:windows
        winupdate=winupdate+1;
        if winupdate>9
        f = waitbar(window/windows,f,'Removing line hum');
        winupdate=0;
        end
    %    fprintf('\b\b\b\b\b%5d',window);
    samples=(1:windowlength)+(window-1)*windowlength;
        fd = d(samples);
    if 0
        if window<2
        [fd,flt]=bandpass(fd,[55 65],Fs); % bandpass around fundamental to make it easier to fit
    else
        fd=filter(flt,fd); % filter already created on first pass
    end
    % fundamental
    end
  [param_fund,fit_fund,fd1,err] = sinefit(fd,time,fundamental,Fs,0,plotting);
%    fd1=removeLineNoise_SpectrumEstimation(fd', Fs);
    % harmonics
    if harmonics
        for h = 1:harmonics
            harm = fundamental*(h+1);
            [~,~,fd1] = sinefit(fd1,time,harm,Fs,0,plotting);
        end
    end
    
    % plot filtered against non-filtered
    if plotting
        figure('name','Raw Vs. Filtered'); 
        plot(time,d,time,fd,'r');
        hold on;
        
        xlabel(['time (',num2str(si),' s)']);
        legend('raw','filtered');
        hold off
    end
    fdout(samples)=fd1; %d(samples)-fit_fund;
    
    end
    close(f);
    fprintf('\n');
    fd=fdout;
end % end of function







