function Project3()
clear
close all
load('data_cn_project_iii_a17.mat')
%
%% 1.Stimulus Nature
corr=xcorr(Stimulus);
t=1/1000:1/1000:20;
figure(1)
plot([-t(50:-1:1) 0 t(1:50)],corr(19950:20050));
xlim([-0.05 0.05])
title('Autocorrelation function of stimulus')
xlabel('\tau')
ylabel('R(\tau)')
%}
%
%% 2.PSTH and mean firing rate
PSTH=zeros(4,20000);
figure(2)
for i=1:4
    for j=1:50
        PSTH(i,:)=PSTH(i,:)+histcounts(All_Spike_Times{i,j}*1000,0:20000)*1000/50;
    end
    subplot(2,2,i)
    stem(t,PSTH(i,:));
    xlabel('time (s)')
    ylabel('Mean Firing Rate (spikes/sec)')
    title(['PSTH - Neuron' num2str(i)])
end

r_pred=PSTH(:,15001:20000);
%}
%
%% 3.Scatter Plots and Poisson Distribution
binsize=[10 20 50 100 200 500];

for k=1:length(binsize)
    spike_cnt=zeros(4,50,20000/binsize(k));
    for i=1:4
        for j=1:50
            spike_cnt(i,j,:)=histcounts(All_Spike_Times{i,j}*1000,0:binsize(k):20000);
        end
        mu=zeros(4,20000/binsize(k));
        sigma2=zeros(4,20000/binsize(k));
        for l=1:20000/binsize(k)
            mu(i,l)=mean(spike_cnt(i,:,l));
            sigma2(i,l)=var(spike_cnt(i,:,l));
        end
        figure(2+k)
        subplot(2,2,i)
        hold on
        scatter(sigma2(i,:),mu(i,:),10)
        plot(mu(i,:),mu(i,:))
        xlabel('\sigma^2')
        ylabel('\mu')
        title(['Neuron ' num2str(i) ' - Binsize ' num2str(binsize(k)) 'ms'])
    end
end
%}
%
%% 4.STA and correction for non-Gaussianity
STA=zeros(4,100);
h=zeros(4,100);
Css=zeros(100,100);

for j=1:14901
        Css=Css+Stimulus(j:j+99)'*Stimulus(j:j+99);
end
Css=Css/14901;
    
for i=1:4
    totalspikes=0;
    for j=1:50
        totalspikes=totalspikes+nnz(All_Spike_Times{i,j}<=15);
        for k=1:nnz(All_Spike_Times{i,j}<=15)
           temp=round(All_Spike_Times{i,j}(k)*1000);
           temp=Stimulus(max([temp-99 1]):temp);
           temp=[zeros(1,100-length(temp)) temp];
           STA(i,:)=STA(i,:)+temp;
        end
    end
    STA(i,:)=STA(i,:)/totalspikes;
    figure(9)
    subplot(2,2,i)
    plot(-99:0,STA(i,:))
    title(['Before Whitening - Neuron ' num2str(i)])
    ylabel('STA')
    xlabel('time (ms)')
    ylim([-0.2 0.2])

    h(i,100:-1:1)=(Css\STA(i,:)')';
    figure(10)
    subplot(2,2,i)
    plot(0:99,h(i,:))
    title(['After Whitening - Neuron ' num2str(i)])
    ylabel('h(t)')
    xlabel('time (ms)')
end
%}
%
%% 5.Determining output non-linearity
y=zeros(4,20099);
binsize=500;
y_binned=zeros(4,15000/binsize);
PSTH_binned=zeros(4,15000/binsize);
figure(11)

for i=1:4
    y(i,:)=conv(Stimulus,h(i,:));
    for j=1:15000/binsize
        y_binned(i,j)=mean(y(i,binsize*(j-1)+1:binsize*j));
        PSTH_binned(i,j)=mean(PSTH(i,binsize*(j-1)+1:binsize*j));
    end
end

ft = fittype( 'a1+a2*x+a3*x^2+a4*x^3', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.186872604554379 0.489764395788231 0.5 0.5];
f=cell(1,4);
for i=1:4
    [xData,yData] = prepareCurveData(y_binned(i,:),PSTH_binned(i,:));
    [f{i},~] = fit(xData,yData,ft,opts);
    subplot(2,2,i)
    plot(f{i},xData,yData);
    legend('hide')
    title(['Non-linearity - Neuron ' num2str(i)])
    xlabel('y(t)')
    ylabel('\lambda(t)')
end
%}
%
%% 6.Prediction Performance and filter pruning
lambda_pred=zeros(4,5000/binsize);
lambda_est=zeros(4,5000/binsize);
figure(12)
for i=1:4
    for j=1:5000/binsize
        x=mean(y(i,15001+binsize*(j-1):15000+binsize*j ));
        lambda_pred(i,j)=f{i}.a1+f{i}.a2*x+f{i}.a3*x.^2+f{i}.a4*x.^3;
        lambda_est(i,j)=mean(r_pred(i,binsize*(j-1)+1:binsize*j));
    end
    subplot(2,2,i)
    scatter(lambda_est(i,:),lambda_pred(i,:),3,'b')
    title(['Prediction - Neuron ' num2str(i)])
    xlabel('\lambda_e_s_t(t)')
    ylabel('\lambda_p_r_e_d(t)')
end

%Filter pruning
figure(13)
h_prune=h;
h_temp=h;
thresh=[0.5 0.48 0.35 0.5];
for i=1:4
    r=corrcoef(lambda_pred(i,:),lambda_est(i,:));
    k=1;
    r2=r(1,2)^2;
    while(k<100)
        k=k+1;
        temp=h_prune(i,:);
        temp(temp==0)=NaN;
        [~,temp]=min(abs(temp));
        h_temp(i,temp)=0;
        y(i,:)=conv(Stimulus,h_temp(i,:));
        for j=1:15000/binsize
            y_binned(i,j)=mean(y(i,binsize*(j-1)+1:binsize*j));
            PSTH_binned(i,j)=mean(PSTH(i,binsize*(j-1)+1:binsize*j));
        end
        [xData,yData] = prepareCurveData(y_binned(i,:),PSTH_binned(i,:));
        [f{i},~]=fit(xData,yData,ft,opts);
        for j=1:5000/binsize
            x=mean(y(i,15001+binsize*(j-1):15000+binsize*j ));
            lambda_pred(i,j)=f{i}.a1+f{i}.a2*x+f{i}.a3*x.^2+f{i}.a4*x.^3;
            lambda_est(i,j)=mean(r_pred(i,binsize*(j-1)+1:binsize*j));
        end
        r=corrcoef(lambda_pred(i,:),lambda_est(i,:));
        r2(k)=r(1,2)^2;
        if(r2(k)<r2(k-1)-thresh(i))
             break
        end
        h_prune(i,temp)=0;
    end
    figure(13)
    subplot(2,2,i)
    plot(100:-1:101-length(r2),r2);
    hold on
    a=plot(102-k,r2(k-1),'+');
    set(a,'linewidth',2)
    hold off
    title(['h Pruning - Neuron ' num2str(i)])
    ylabel('r^2')
    xlabel('No. of parameters in h(t)')
end

for i=1:4
    figure(14)
    subplot(2,2,i)
    plot(0:99,h_prune(i,:))
    title(['Pruned filter - Neuron ' num2str(i)])
    xlabel('time (ms)')
    ylabel('h(t)')
    figure(15)
    subplot(2,2,i)
    plot(abs(fft(h(i,:))))
    hold on
    plot(abs(fft(h_prune(i,:))))
    hold off
    title(['FFT of h(t) - Neuron ' num2str(i)])
    legend('Normal filter','Pruned filter')
end
figure(13)
%}
%
%% 7.B.Discrimination based on VP-SDM
q=[0 0.001 0.01 0.1 1 10 100];
MI=zeros(4,length(q));
for iter=1:100
    while(1)
    	idx=randsample(19901,8)';               %random indices
        a=abs(bsxfun(@minus,idx,idx'));
        if(min(a(a>0))>=100)                    %check if non-overlapping
            break
        end
    end

    for m=1:4
        t=cell(8,50);
        for i=1:8
            for j=1:50
                temp=All_Spike_Times{m,j};
                t{i,j}=temp(temp>=idx(i)/1000&temp<(idx(i)+100)/1000);
            end
        end

        for n=1:length(q)
            confusion=zeros(8,8);
            for i=1:8
                for j=1:50
                    mean_dist=zeros(1,8);
                    for k=1:8
                        for l=1:50
                            if(k==i&&l==j)
                                continue;
                            end
                            mean_dist(k)=mean_dist(k)+VPSDM(t{i,j},t{k,l},q(n));
                        end
                    end
                    mean_dist=mean_dist/50;
                    mean_dist(i)=mean_dist(i)*50/49;
                    [~,k]=min(mean_dist);
                    confusion(i,k)=confusion(i,k)+1;
                end
            end
            confusion=confusion/50;
            MI(m,n)=MI(m,n)+MutualInfo(confusion)/100;
        end
    end
end

figure(16)
for m=1:4
    subplot(2,2,m)
    errorbar(-3:2,MI(m,2:end),0.05*MI(m,2:end))
    hold on
    [~,p]=max(MI(m,:));
    p=plot(log10(q(p)),MI(m,p),'+');
    set(p,'linewidth',2)
    title(['Discrimination - Neuron ' num2str(m)])
    xlabel('log_1_0(q)')
    ylabel('MI(q)')
end
%}
end

function d=VPSDM(tli,tlj,q)
nspi=length(tli);
nspj=length(tlj);

if q==0
   d=abs(nspi-nspj);
   return
elseif q==Inf
   d=nspi+nspj;
   return
end

scr=zeros(nspi+1,nspj+1);
scr(:,1)=(0:nspi)';
scr(1,:)=(0:nspj);
if(nspi && nspj)
   for i=2:nspi+1
      for j=2:nspj+1
         scr(i,j)=min([scr(i-1,j)+1 scr(i,j-1)+1 scr(i-1,j-1)+q*abs(tli(i-1)-tlj(j-1))]);
      end
   end
end
d=scr(nspi+1,nspj+1);
end

function MI = MutualInfo(confusion)
MI=0;
for i=1:size(confusion,1)
    for j=1:size(confusion,2)
        if(confusion(i,j)~=0)
            MI=MI+confusion(i,j)/size(confusion,1)*log2(confusion(i,j)/sum(confusion(:,j)));          %confusion matrix has entries of p(y/x)
        end
    end
end
end