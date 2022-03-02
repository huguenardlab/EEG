function openephysreader(filePath,varargin)
%
%                    >>> INPUT VARIABLES >>>
%
% NAME        TYPE, DEFAULT      DESCRIPTION
% fn          char array         path where abf files to be analyzed are located
% All optional input parameters listed below (= all except the file name)
% must be specified as parameter/value pairs, e.g. as in
%          openephysreader('/home/shared/Astra/EEG','lowpassfilter',25,'highpassfilter',1);

inputpath=filePath;
outputpath=inputpath;
framewidth=5; % in seconds
EEGmax=.4
use_previous=1; % use previous manually detected events if available.
show_wavelet_events=1;  % display seizures auto detected
show_rejected_events=1;
% looking for abf files.  If there are some ABF files, then the path1 statement below will need to be changed.
exportfmt='-depsc2';

clipping=.6; % mV level to get rid of movement artifacts, etc
artifact_suppression=0; % number of standard deviations from mean to truncate the data

autoquit=0;


% the following are values to use for a straight xltek->abf converted file
% this will bandpass filter between 1 and 25 Hz, and resample the original 500Hz sampling
% rate down to 100 Hz, which is plenty for looking at SWD
% frequency in Hz
transsuppression=0; %mv in sudden change in eeg trace to consider it an artifact
maxtrans=1000; %maximum duration of such an artifact to consider eliminating it
medianfilter=0;

absolutethreshold=0;
lowpassfilter=0;
highpassfilter=1;
filtering=1;
MaxSampleRate=5000; %downsample to this frequency before proceeding

saving_episodes=1;
pausing=0;
previewing=0;
if previewing>0
    usedchannels=[4];
    pausing=-1;
    num_plotted_events=1e6;
    saving_episodes=0;
end
printing_summary=1;
fs=250; % 250 Hz should be plenty for EEG, but can over-ride if desired with input arguments


%END OF VARIABLES
% assign values of optional input parameters, if any were given
pvpmod (varargin);
numfiles=1;
disp(filePath);
d = dir (filePath);
numfiles=length(d);
if numfiles<1
    disp('No files found');
end
fnames=[];
nfolders=0;
for i=1:numfiles
    if d(i).isdir==1 && ~strcmp(d(i).name,'.') && ~strcmp(d(i).name,'..')
        nfolders=nfolders+1;
        fnames{nfolders}=d(i).name;
    end
end
figure (1);
plot(0,0);
datacursormode off;
numfiles=size(fnames,2);
for thisfile = 1:numfiles
    tic;
    fna=fnames{thisfile}; 
    % the following code will decode file names to get cohort, subject, time
    % and date
    trialid=split(fna,'_');
    trialsz=size(trialid,1);
    cohort=trialid{1};
    tdate=trialid{trialsz-1};
    ttime=trialid{trialsz};
    nsubjects=1;
    subid=split(trialid{2},'.');
    subn=1;
    subj(subn)=str2num(subid{1});
    chs(subn,:)=[str2num(subid{2}) str2num(subid{3})];
    if (trialsz>4)
        % this section is hard coded to deal with file name format that might
        % have two animals in an individual EEG file
        subn=2;
        subid=split(trialid{3},'.');
        subj(subn)=str2num(subid{1});
        chs(subn,:)=[str2num(subid{2}) str2num(subid{3})];
        nsubjects=2;
    end
    num_kept_channels=2;kept_channels=chs(subn,:);
    thisfile=[filePath '/' fna];
    s=sprintf('Reading File: %s, Cohort: %s, Date %s, Time %s',fna,cohort,tdate,ttime);
    disp(s);
    [fileData,fileHeader] = OpenEphysLoader( thisfile,'desiredfs',fs); % load EEG channels
    
    for subject=1:nsubjects
        numsz=0;
        szstart=[];
        szend=[];
        
        d=fileData;
        si=1000/fs; %ms;
        
        if filtering>0
            if lowpassfilter > 0
                [blpf,alpf]=butter(2,lowpassfilter*2*si/1000);
            end
            if highpassfilter >0
                [bhpf,ahpf]=butter(2,highpassfilter*2*si/1000,'high');
            end
        end
        s=sprintf('Analyzing Subject %d',subj(subject));
        disp(s);
        chans=[];
        for thisch=1:size(chs,2)
            chans=[chans find(fileHeader.chans==chs(subject,thisch),1)]; % find first occurence
        end
        sst=d(:,chans)/1000; % divide by 1000 to put in mV 

        if filtering>0
            % band pass filter
            if lowpassfilter>0
                sst=filter(blpf,alpf,sst);
            end
            if highpassfilter>0
                sst=filter(bhpf,ahpf,sst);
            end
        end
        % or even 50 Hz does a fine job, but must adjusct s0 accordingly below
        time=0:length(sst)-1;
        time=time.*si;
        times=time/1000.;
        n = length(sst);
        dt=time(2)-time(1);
        done=0;
        
        frame=0;
        framepoints=framewidth/si*1000;
        insz=0;
        %------------------------------------------------------ Plotting
        %--- Plot time series
        szstartar=[];
        quitting=0;
        ttl=sprintf('File %s, Subject %d',fna,subj(subject));
        fout=sprintf('%s/%s_subject_%d_seizures.txt',outputpath,fna,subj(subject));
        if show_wavelet_events
            fout1=sprintf('%s/%s_subject%d_ch%d_results.txt',outputpath,fna,subj(subject),chs(subject,1));
            fout2=sprintf('%s/%s_subject%d_ch%d_results.txt',outputpath,fna,subj(subject),chs(subject,2));
            fnr1=fopen(fout1,'r');
            fnr2=fopen(fout2,'r');
            if fnr1    >0% outputfile exists
                fgetl(fnr1); % get and toss first header line
                c=textscan(fnr1,'%f');
                tmp=c{1};
                nrows=size(tmp,1)/8;
                ncols=size(tmp,1)/nrows;
                tmp2=reshape(tmp,ncols,nrows)';
                szstartar(1:nrows,1)=tmp2(:,1);
                szendar(1:nrows,1)=tmp2(:,2);
                szrejar(1:nrows,1)=tmp2(:,8);
                fclose(fnr1);
            end
            if fnr2    >0% outputfile exists
                fgetl(fnr2); % get and toss first header line
                c=textscan(fnr2,'%f');
                tmp=c{1};
                nrows=size(tmp,1)/8;
                ncols=size(tmp,1)/nrows;
                tmp2=reshape(tmp,ncols,nrows)';
                szstartar(1:nrows,2)=tmp2(:,1);
                szendar(1:nrows,2)=tmp2(:,2);
                szrejar(1:nrows,2)=tmp2(:,8);
                fclose(fnr2);
            end
        end
        if use_previous
            fnr=fopen(fout,'r');
            if fnr    >0% outputfile exists
                fgetl(fnr); % get and toss first two header lines
                fgetl(fnr);
                c=textscan(fnr,'%f %f');
                szstart=c{1};
                szend=c{2};
                fclose(fnr);
            end
            
        end
        prompt='''N''ext, ''P''revious, ''J''ump, ''Q''uit';
        nbytes=fprintf('%s',prompt);
        while ~done
            firstpt=frame*framepoints+1;
            lastpt=firstpt+framepoints-1;
            
            plot (times(firstpt:lastpt),sst(firstpt:lastpt,1),'k');
            title(ttl,'Interpreter','none')
            xlabel('Time (s)');
            ylabel('EEG (mV)');
            hold on;
            plot (times(firstpt:lastpt),sst(firstpt:lastpt,2)-2*EEGmax,'k');
            ylim([-EEGmax*3,EEGmax]);
            xlim([times(firstpt) times(lastpt)]);
            thisframedone=0;
            
            while ~thisframedone
                thesestarts= (szstart(szstart>=times(firstpt) &szstart<=times(lastpt)));
                for s=1:length(thesestarts)
                    plot([thesestarts(s),thesestarts(s)],[-3*EEGmax,EEGmax],'r','LineWidth',2);
                    text(thesestarts(s),-EEGmax,'\leftarrow Manual seizure start');
                end
                theseends= (szend(szend>=times(firstpt) &szend<=times(lastpt)));
                for s=1:length(theseends)
                    plot([theseends(s),theseends(s)],[-3*EEGmax,EEGmax],'b','LineWidth',2);
                    text(theseends(s),-EEGmax,'\leftarrow Manual seizure end');
                end
                if show_wavelet_events && size(szstartar,1)>0
                    for thisch=1:2
                        szstarta=nonzeros(szstartar(:,thisch));
                        szenda=nonzeros(szendar(:,thisch));
                        szreja=szrejar(:,thisch);                        
                    thesestartsa= (szstarta(szstarta>=times(firstpt) &szstarta<=times(lastpt)));
                for s=1:length(thesestartsa)
                    idx=find(szstarta==thesestartsa(s));
                    
                        rej1='';
                        if szreja(idx)>0 
                            rej1=' REJECTED';
                            rejcolor=[.6 0 0];
                            not_rejected=0;
                        else
                            rejcolor='g';
                            not_rejected=1;
                        end
                if not_rejected || show_rejected_events                        
                    plot([thesestartsa(s),thesestartsa(s)],[-EEGmax,EEGmax]-(thisch-1)*EEGmax*2,'color',rejcolor,'LineWidth',not_rejected+1);
                    text(thesestartsa(s),-EEGmax-(thisch-1.5)*EEGmax*1.5,['\leftarrow DETECTED seizure start' rej1],'color',rejcolor);
                    end
                end
                theseendsa= (szenda(szenda>=times(firstpt) &szenda<=times(lastpt)));
                for s=1:length(theseendsa)
                     idx=find(szenda==theseendsa(s));
                        rej1='';
                        if szreja(idx)>0 
                            rej1=' REJECTED';
                            rejcolor='c';
                            not_rejected=0;
                        else
                            not_rejected=1;
                            rejcolor='m';
                        end
            
                    if  not_rejected || show_rejected_events                    
                        plot([theseendsa(s),theseendsa(s)],[-EEGmax,EEGmax]-(thisch-1)*EEGmax*2,'color',rejcolor,'LineWidth',not_rejected+1);
                        text(theseendsa(s),-EEGmax-(thisch-1.5)*EEGmax*.8,['\leftarrow DETECTED seizure end' rej1],'color',rejcolor);
                    end
                end
                    end
                    
                end
                prompt='''N''ext, ''P''revious, ''J''ump, ''W''iden, ''T''ighten, ''U''p, ''D''own, ''Q''uit';
                text((times(firstpt)+times(lastpt))/2, EEGmax*1.25,prompt,'HorizontalAlignment','center');
                if insz
                    prompt1=['''E''nd of seizure,   ' prompt];
            
                else
                    
                    prompt1=['''S''tart of seizure, ' prompt];
                end
                
                fprintf(repmat('\b',1,nbytes));
                nbytes=fprintf('%s',prompt1);
                k=waitforbuttonpress;
                if k
                    cchr=get(gcf,'CurrentCharacter');
                else
                    cchr='{'; % mouse clicked.   don't return character.
                end
                if cchr=='q'
                    quitting=1;
                    thisframedone=1;
                end
                if cchr=='s' && ~insz
                    fprintf(repmat('\b',1,nbytes));
                    nbytes=fprintf('Click on screen to mark onset time\n');
                    [x,y]=ginput(1);
                    numsz=numsz+1;
                    insz=1;
                    szstart(numsz)=x;
                end
                if cchr=='e' && insz
                    fprintf(repmat('\b',1,nbytes));
                    nbytes=fprintf('Click on screen to mark offset time\n');
                    [x,y]=ginput(1);
                    szend(numsz)=x;
                    insz=0;
                end
                if cchr=='j'
                    fprintf(repmat('\b',1,nbytes));                    
                    timejump=input('Time to Jump to (s): ','s');
                    nbytes=length('Time to Jump to (s): ')+length(timejump)+1;
                    frame=floor(str2double(timejump)/framewidth);
                    thisframedone=1;
                end
                if cchr=='n'
                    frame=frame+1;
                    thisframedone=1;
                end
                if cchr=='p'
                    if frame>0;
                        frame=frame-1;
                    end
                    thisframedone=1;
                end
                if cchr=='w'
                    framewidth=framewidth*2;
                    framepoints=framewidth/si*1000;
                    frame=frame/2;
                    thisframedone=1;
                end
                if cchr=='t'
                    framewidth=(framewidth/2);
                    framepoints=floor(framewidth/si*1000);
                    frame=frame*2;
                    thisframedone=1;                    
                end
                if cchr=='u'
                    EEGmax=EEGmax/2;
                    thisframedone=1;                    
                end
                if cchr=='d'
                    EEGmax=EEGmax*2;                    
                    thisframedone=1;
                end

                if quitting || frame*framepoints+framepoints>size(sst,1)
                    done=1;
                end
            end
            hold off;
        end
        
        fnr=fopen(fout,'w');
        fprintf(fnr,'seizuresevents from file %s.abf subject %d\nstart\tszend\n',fna,subj(subject));
        for s=1:numsz
            fprintf(fnr,'%f\t%f\n',szstart(s),szend(s));
        end
        fclose(fnr);
        
    end
    
    TimeN=toc;
    fprintf('Time = %3.3f sec\n', TimeN);
end %fname file list of abf files
fprintf('Done!\n');

if autoquit >0
    quit;
end
