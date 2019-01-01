close all;
clear;

%%%%---- THIS MARKS THE START OF PART A OF THE PROJECT ----%%%%

%{
%% 1. Responses to Tones (rate representation)
cohc  = 1.0;   % normal ohc function
cihc  = 1.0;   % normal ihc function
fiberType = 3; % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
implnt = 0;    % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse

Fs = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 50e-3;  % stimulus duration in seconds
rt = 5e-3;   % rise/fall time in seconds

% PSTH parameters
nrep = 1;

t = 0:1/Fs:T-1/Fs; % time vector
mxpts = length(t);
irpts = rt*Fs;

BF=[.5e3 4e3];
f=125*2.^(0:1/8:8);
tuning.BF500=zeros(10,length(f));
tuning.BF4000=zeros(10,length(f));
stimdb = -10:10:80;

for k=1:length(stimdb)
    for i=1:length(f)
        pin = sqrt(2)*20e-6*10^(stimdb(k)/20)*sin(2*pi*f(i)*t);    % Single tone at frequency f(i)
        pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;
        pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
        for j=1:20
            vihc = catmodel_IHC(pin,BF(1),nrep,1/Fs,T*2,cohc,cihc);
            [synout,~] = catmodel_Synapse(vihc,BF(1),nrep,1/Fs,fiberType,implnt);
            tuning.BF500(k,i)=tuning.BF500(k,i)+max(synout)/20;
            vihc = catmodel_IHC(pin,BF(2),nrep,1/Fs,T*2,cohc,cihc);
            [synout,~] = catmodel_Synapse(vihc,BF(2),nrep,1/Fs,fiberType,implnt);
            tuning.BF4000(k,i)=tuning.BF4000(k,i)+max(synout)/20;
        end
    end
end

figure(1)
hold on
for k=1:length(stimdb)
    plot(log2(f),tuning.BF500(k,:))
end
xlim([log2(125) log2(32000)]);
set(gca,'XTick',log2(125*2.^(0:8)),'XTickLabel',125*2.^(0:8));
title(['ANF tuning curve for BF = ' num2str(BF(1)) ' Hz']);
xlabel('f (Hz)');
ylabel('Firing rate (spikes/sec)')
legend('-10dB SPL','0dB SPL','10dB SPL','20dB SPL','30dB SPL','40dB SPL','50dB SPL','60dB SPL','70dB SPL','80dB SPL');
figure(2)
hold on
for k=1:length(stimdb)
    plot(log2(f),tuning.BF4000(k,:))
end
xlim([log2(125) log2(32000)]);
set(gca,'XTick',log2(125*2.^(0:8)),'XTickLabel',125*2.^(0:8));
title(['ANF tuning curve for BF = ' num2str(BF(2)) ' Hz']);
xlabel('f (Hz)');
ylabel('Firing rate (spikes/sec)')
legend('-10dB SPL','0dB SPL','10dB SPL','20dB SPL','30dB SPL','40dB SPL','50dB SPL','60dB SPL','70dB SPL','80dB SPL');

ratevsint=zeros(2,length(stimdb));
ratevsint(1,:)=tuning.BF500(:,find(f==BF(1)));
ratevsint(2,:)=tuning.BF4000(:,find(f==BF(2)));
figure(3)
hold on
for i=1:2
	plot(stimdb,ratevsint(i,:))
end
xlabel('Intensity (dB SPL)')
ylabel('Firing rate at BF (spikes/sec)')
title('Rate vs. Intensity plot')
legend('BF = 500Hz','BF = 4kHz')
%}
%{
%% 2. Responses to Speech (rate representation)
% Steady state calibration
BF=125*2.^(0:1/8:6);
[speech.orig,~]=audioread('fivewo.wav');
audiowrite('speech.wav',speech.orig,100000);
[speech.orig,Fs]=audioread('speech.wav');
speech.steady=speech.orig(105001:115000);
speech.int=20*log10(rms(speech.steady)/20e-6);

stimdb = -20:10:80;
nrep = 1;
mxpts = length(speech.steady);
T=length(speech.steady)/Fs;
ratevsint=zeros(1,length(stimdb));

for i=1:length(stimdb)
    pin = speech.steady'*10^((stimdb(i)-speech.int)/20);
    pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;
    pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
    for j=1:20
        vihc = catmodel_IHC(pin,500,nrep,1/Fs,T,cohc,cihc);         % No resting period taken
        [synout,~] = catmodel_Synapse(vihc,500,nrep,1/Fs,fiberType,implnt);
        ratevsint(1,i)=ratevsint(1,i)+max(synout)/20;
    end
end

figure(4)
plot(stimdb,ratevsint)
xlabel('Intensity (dB SPL)')
ylabel('Firing rate at BF (spikes/sec)')
title('Rate vs. Intensity plot for BF = 500Hz')

% Spiking response to entire speech signal
stimdb=[20 50 80];              % 3 sound levels seen from plot
T=length(speech.orig)/Fs;
mxpts = length(speech.orig);
nrep=50;
spiketrain=cell(length(stimdb),length(BF));

for i=1:length(stimdb)
    pin = speech.orig'*10^((stimdb(i)-speech.int)/20);
    pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;
    pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
    for j=1:length(BF)
        vihc = catmodel_IHC(pin,BF(j),nrep,1/Fs,T*1.1,cohc,cihc);
        [~,spiketrain{i,j}] = catmodel_Synapse(vihc,BF(j),nrep,1/Fs,fiberType,implnt);
        spiketrain{i,j}(length(speech.orig)+1:end)=[];
    end
end

% Spectrogram analysis
binlen=[4e-3 8e-3 16e-3 32e-3 64e-3 128e-3]*Fs;
for l=1:length(stimdb)
    figure(3+2*l)
    winlen=25.6e-3*Fs;
    spectrogram(speech.orig*10^((stimdb(l)-speech.int)/20),hann(winlen),winlen/2,winlen,Fs,'yaxis');
    title(['Spectrogram of Speech signal (' num2str(stimdb(l)) 'dB SPL)']);

    figure(4+2*l)
    for i=1:length(binlen)
        spect=zeros(length(BF),floor(2*length(spiketrain{l,1})/binlen(i))-1);
        for j=1:length(BF)
            for k=1:size(spect,2)
                spect(j,k)=sum(spiketrain{l,j}((k-1)*binlen(i)/2+1:(k+1)*binlen(i)/2))*Fs/(nrep*binlen(i)); %moving window with 50% overlap
            end
        end
        subplot(3,2,i)
        image([0 T],[log2(BF(1)) log2(BF(end))],flipud(spect),'CDataMapping','scaled');
        set(gca,'YTick',log2(125*2.^(0:6)),'YTickLabel',125*2.^(6:-1:0)/1000);
        xlabel('Time (secs)')
        ylabel('BF of ANFs (kHz)')
        title(['Rate proxy spectrogram: ' num2str(stimdb(l)) 'dB SPL w/ window = ' num2str(binlen(i)*1000/Fs) 'ms'])
        colorbar;
    end
end
%}
%{
%% 3. Responses to Speech (fine timescale representation)
for l=1:3
    binlen=0.1e-3*Fs;
    BF=125*2.^(0:1/2:5);
    psth=zeros(length(BF),floor(2*length(spiketrain{l,1})/binlen)-1);
    for j=1:length(BF)
        for k=1:size(psth,2)
            psth(j,k)=sum(spiketrain{l,j}((k-1)*binlen/2+1:(k+1)*binlen/2))*Fs/(nrep*binlen);
        end
    end

    winlen=12.8e-3*Fs;
    domf=zeros(length(BF),floor(2*size(psth,2)/winlen)-1);
    figure(10+l)
    subplot(2,1,1)
    spectrogram(speech.orig*10^((stimdb(l)-speech.int)/20),hann(binlen),binlen/2,binlen,Fs,'yaxis');
    title(['Spectrogram of Speech signal (' num2str(stimdb(l)) 'dB SPL)']);
    ylim([0 4])

    t=(winlen/2:winlen/2:size(psth,2)-winlen/2);
    t=(t+1)*binlen/2/Fs;                            % getting t in secs
    
    subplot(2,1,2)
    ylim([0 4])
    hold on
    for j=1:length(BF)
        for k=1:size(domf,2)
            temp=abs(fft(psth(j,(k-1)*winlen/2+1:(k+1)*winlen/2)));
            temp=temp(2:end/2+1);
            [~,domf(j,k)]=max(temp);
            domf(j,k)=domf(j,k)/length(temp)*Fs;      % dominant frequency
        end
        plot(t,domf(j,:)/1000,'*')
    end
    title('Dominant frequencies')
    h=get(gca,'Children');
    legend([h(11) h(10) h(9) h(8) h(7) h(6) h(5) h(4) h(3) h(2) h(1)],['BF=' num2str(BF(1)) 'Hz'],['BF=' num2str(BF(2)) 'Hz'],['BF=' num2str(BF(3)) 'Hz'],['BF=' num2str(BF(4)) 'Hz'],['BF=' num2str(BF(5)) 'Hz'],['BF=' num2str(BF(6)) 'Hz'],['BF=' num2str(BF(7)) 'Hz'],['BF=' num2str(BF(8)) 'Hz'],['BF=' num2str(BF(9)) 'Hz'],['BF=' num2str(BF(10)) 'Hz'],['BF=' num2str(BF(11)) 'Hz'])
end
%}
%%
%%%%---- THIS MARKS THE START OF PART B OF THE PROJECT ----%%%%
%{
%% Reconstructing sound using primarily temporal cues
close all;
nbands=[1 2 8 4];
speech.reconstruct=zeros(4,length(speech.orig));
for j=1:length(nbands)
    filtcent=250*2.^(0:3/(nbands(j)-1):3);
    speech.filt=zeros(nbands(j),length(speech.orig));
    fc=0.004/nbands(j);
    
    for i=1:nbands(j)
        [b,a]=butter(2,[filtcent(i)*2/Fs-fc filtcent(i)*2/Fs+fc],'bandpass');
        speech.filt(i,:)=filtfilt(b,a,speech.orig);
        speech.filt(i,:)=abs(hilbert(speech.filt(i,:)));    % Hilbert transform
        [b,a]=butter(4,50*2/Fs,'low');
        speech.filt(i,:)=filtfilt(b,a,speech.filt(i,:));    % Obtaining the envelope
        speech.filt(i,:)=speech.filt(i,:).*wgn(1,size(speech.filt,2),50); % Modulating white gaussian noise with envelope
        [b,a]=butter(2,[filtcent(i)*2/Fs-fc filtcent(i)*2/Fs+fc],'bandpass');
        speech.filt(i,:)=filtfilt(b,a,speech.filt(i,:));  % Limiting frequencies of modulated noise
    end

    speech.reconstruct(j,:)=sum(speech.filt,1);
    [b,a]=butter(4,4e3*2/Fs,'low');
    speech.reconstruct(j,:)=filtfilt(b,a,speech.reconstruct(j,:));
end
%}
%{
%% Repetition of part A2
% Steady state calibration
BF=125*2.^(0:1/8:6);
for m=1:3:4
    speech.steady=speech.reconstruct(m,105001:115000);
    speech.int=20*log10(rms(speech.steady)/20e-6);

    stimdb = -20:10:80;
    nrep = 1;
    mxpts = length(speech.steady);
    T=length(speech.steady)/Fs;
    ratevsint=zeros(1,length(stimdb));

    for i=1:length(stimdb)
        pin = speech.steady*10^((stimdb(i)-speech.int)/20);
        pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;
        pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
        for j=1:20
            vihc = catmodel_IHC(pin,500,nrep,1/Fs,T,cohc,cihc);         % No resting period taken
            [synout,~] = catmodel_Synapse(vihc,500,nrep,1/Fs,fiberType,implnt);
            ratevsint(1,i)=ratevsint(1,i)+max(synout)/20;
        end
    end
    if(m==1)
        figure(14)
    elseif(m==4)
        figure(21)
    end
    plot(stimdb,ratevsint)
    xlabel('Intensity (dB SPL)')
    ylabel('Firing rate at BF (spikes/sec)')
    title(['Rate vs. Intensity plot: BF = 500Hz, (' num2str(m) ' bands)'])

    % Spiking response to entire speech signal
    stimdb=[20 50 80];              % 3 sound levels seen from plot
    T=size(speech.reconstruct,2)/Fs;
    mxpts = size(speech.reconstruct,2);
    nrep=50;
    spiketrain=cell(length(stimdb),length(BF));

    for i=1:length(stimdb)
        pin = speech.reconstruct(m,:)*10^((stimdb(i)-speech.int)/20);
        pin(1:irpts)=pin(1:irpts).*(0:(irpts-1))/irpts;
        pin((mxpts-irpts):mxpts)=pin((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;
        for j=1:length(BF)
            vihc = catmodel_IHC(pin,BF(j),nrep,1/Fs,T*1.1,cohc,cihc);
            [~,spiketrain{i,j}] = catmodel_Synapse(vihc,BF(j),nrep,1/Fs,fiberType,implnt);
            spiketrain{i,j}(size(speech.reconstruct,2)+1:end)=[];
        end
    end

    % Spectrogram analysis
    binlen=[4e-3 8e-3 16e-3 32e-3 64e-3 128e-3]*Fs;
    for l=1:length(stimdb)
        if(m==1)
            figure(13+2*l)
        elseif(m==4)
            figure(20+2*l)
        end
        winlen=25.6e-3*Fs;
        spectrogram(speech.reconstruct(m,:)*10^((stimdb(l)-speech.int)/20),hann(winlen),winlen/2,winlen,Fs,'yaxis');
        title(['Spectrogram of Speech signal (' num2str(m) ' bands) (' num2str(stimdb(l)) 'dB SPL)']);
        ylim([0 4])
        
        if(m==1)
            figure(14+2*l)
        elseif(m==4)
            figure(21+2*l)
        end
        suptitle([ num2str(m) ' frequency bands'])
        
        for i=1:length(binlen)
            spect=zeros(length(BF),floor(2*length(spiketrain{l,1})/binlen(i))-1);
            for j=1:length(BF)
                for k=1:size(spect,2)
                    spect(j,k)=sum(spiketrain{l,j}((k-1)*binlen(i)/2+1:(k+1)*binlen(i)/2))*Fs/(nrep*binlen(i)); %moving window with 50% overlap
                end
            end
            subplot(3,2,i)
            image([0 T],[log2(BF(1)) log2(BF(end))],flipud(spect),'CDataMapping','scaled');
            set(gca,'YTick',log2(125*2.^(0:6)),'YTickLabel',125*2.^(6:-1:0)/1000);
            xlabel('Time (secs)')
            ylabel('BF of ANFs (kHz)')
            title(['Rate proxy spectrogram: ' num2str(stimdb(l)) 'dB SPL w/ window = ' num2str(binlen(i)*1000/Fs) 'ms'])
            colorbar;
        end
    end
end
%}
%{
%% 3. Repetition of part A3
for m=1:3:4
    for l=1:3
        binlen=0.1e-3*Fs;
        BF=125*2.^(0:1/2:5);
        psth=zeros(length(BF),floor(2*length(spiketrain{l,1})/binlen)-1);
        for j=1:length(BF)
            for k=1:size(psth,2)
                psth(j,k)=sum(spiketrain{l,j}((k-1)*binlen/2+1:(k+1)*binlen/2))*Fs/(nrep*binlen);
            end
        end

        winlen=12.8e-3*Fs;
        domf=zeros(length(BF),floor(2*size(psth,2)/winlen)-1);
        if(m==1)
            figure(27+l)
        elseif(m==4)
            figure(30+l)
        end
        subplot(2,1,1)
        spectrogram(speech.reconstruct(m,:)*10^((stimdb(l)-speech.int)/20),hann(binlen),binlen/2,binlen,Fs,'yaxis');
        title(['Spectrogram of Speech signal (' num2str(m) ' bands, ' num2str(stimdb(l)) 'dB SPL)']);
        ylim([0 4])

        t=(winlen/2:winlen/2:size(psth,2)-winlen/2);
        t=(t+1)*binlen/2/Fs;                            % getting t in secs

        subplot(2,1,2)
        ylim([0 4])
        hold on
        for j=1:length(BF)
            for k=1:size(domf,2)
                temp=abs(fft(psth(j,(k-1)*winlen/2+1:(k+1)*winlen/2)));
                temp=temp(2:end/2+1);
                [~,domf(j,k)]=max(temp);
                domf(j,k)=domf(j,k)/length(temp)*Fs;      % dominant frequency
            end
            plot(t,domf(j,:)/1000,'*')
        end
        title('Dominant frequencies')
        xlabel('Time (secs)')
        ylabel('Frequency (kHz)')
        h=get(gca,'Children');
        legend([h(11) h(10) h(9) h(8) h(7) h(6) h(5) h(4) h(3) h(2) h(1)],['BF=' num2str(BF(1)) 'Hz'],['BF=' num2str(BF(2)) 'Hz'],['BF=' num2str(BF(3)) 'Hz'],['BF=' num2str(BF(4)) 'Hz'],['BF=' num2str(BF(5)) 'Hz'],['BF=' num2str(BF(6)) 'Hz'],['BF=' num2str(BF(7)) 'Hz'],['BF=' num2str(BF(8)) 'Hz'],['BF=' num2str(BF(9)) 'Hz'],['BF=' num2str(BF(10)) 'Hz'],['BF=' num2str(BF(11)) 'Hz'])
    end
end
%}