function main()
clear
close all
%% 1.Making the VerNOP model
% Inhomogenous Poisson
t=1:1000;
lambda=45+2*randn(1,1000);
n=640;
spikes=zeros(n,1000);
for i=1:n
	spikes(i,:)=binornd(ones(1,length(lambda)),lambda/1000);
end
figure(1)
subplot(2,1,1)
plot(t/1000,lambda,'b')
title('Inhomogenous Poisson Process')
xlabel('t(s)')
ylabel('\lambda(t)')
subplot(2,1,2)
plot(t/1000,sum(spikes(1:100,:))*1000/100,'r')
title('Generated PSTH')
xlabel('t(s)')
ylabel('Mean firing rate (spikes/sec)')

n=[10 20 40 80 160 320 640];
rmse=zeros(1,length(n));
figure(2)
for i=1:length(n)
    rmse(i)=sqrt(mean((sum(spikes(1:n(i),:))*1000/n(i)-lambda).^2));
end
plot(n,rmse,'bo-');
title('Error in PSTH')
ylabel('RMSE')
xlabel('n')
lambda=[];

% Stimulus Generation
stream=char(['S'*ones(1,7) 'D' 'S'*ones(1,7)]);
lambda.S=0.5*ones(1,300*length(stream));
lambda.D=0.5*ones(1,300*length(stream));
for i=1:length(stream)
   if(stream(i)=='S')
       lambda.S(300*(i-1)+1:300*(i-1)+50)=10;
       lambda.D(300*(i-1)+1:300*(i-1)+50)=2.5;
   elseif(stream(i)=='D')
       lambda.S(300*(i-1)+1:300*(i-1)+50)=2.5;
       lambda.D(300*(i-1)+1:300*(i-1)+50)=10;        
   end
end

% Model parameters
w.ThS.SP=0.2;
w.ThD.SP=0.2;
w.ThS.L4=0.02;
w.ThD.L4=0.02;
w.SP.L4=0.11;
tau.re=0.9;
tau.ei.Th=10;
tau.ei.SP=27;
tau.ir=5000;
tau.syn=10;
tau.ref=2;

spike.ThS=binornd(ones(1,length(lambda.S)),lambda.S/1000);
spike.ThD=binornd(ones(1,length(lambda.D)),lambda.D/1000);
spike1=spike;                           % Used later for comparison

[V,spike]=verNOP(spike,tau,w);      % No Long-term plasticity

figure(3)
t=1e-3:1e-3:4.5;
subplot(4,1,1)
plot(t,spike.ThS)
title('1. Spiking activity of neurons')
ylabel('Spike.ThS')
subplot(4,1,2)
plot(t,spike.ThD)
ylabel('Spike.ThD')
subplot(4,1,3)
plot(t,spike.SP)
ylabel('Spike.SP')
subplot(4,1,4)
plot(t,spike.L4)
ylabel('Spike.L4')
xlabel('t(s)')
for i=1:4
   subplot(4,1,i)
   ylim([0 1])
end
figure(4)
subplot(2,1,1)
plot(t,V.sp)
title('1. Neuron membrane potential')
ylabel('V_S_P(t)')
subplot(2,1,2)
plot(t,V.l4)
ylabel('V_L_4(t)')
xlabel('t(s)')

%% 2,3.Plotting PSTH and seeing effect of changing tau_ir
temp=[1000 3000 10000 5000];
for j=1:4
PSTH.SP=zeros(1,length(lambda.S)/10);
PSTH.L4=zeros(1,length(lambda.S)/10);
tau.ir=temp(j);
for i=1:50
    spike.ThS=binornd(ones(1,length(lambda.S)),lambda.S/1000);
    spike.ThD=binornd(ones(1,length(lambda.D)),lambda.D/1000);
    [~,spike]=verNOP(spike,tau,w);
    spike_times=find(spike.SP); 
    PSTH.SP=PSTH.SP+histcounts(spike_times,0:10:length(lambda.S))*1000/50/10;
    spike_times=find(spike.L4);
    PSTH.L4=PSTH.L4+histcounts(spike_times,0:10:length(lambda.S))*1000/50/10;
end

figure(4+j)
subplot(2,1,1)
stem(PSTH.SP)
title(['PSTH of SP Neuron (\tau_i_r=' num2str(temp(j)) 'ms)'])
subplot(2,1,2)
stem(PSTH.L4)
title(['PSTH of L4 Neuron (\tau_i_r=' num2str(temp(j)) 'ms)'])
for i=1:2
    subplot(2,1,i)
    ylabel('PSTH(sp/s)')
end
end

%% 4.a.Long Term Plasticity
stream=char('S'*ones(1,12000));
stream(binornd(ones(1,12000),0.1*ones(1,12000))==1)='D'; %Probability of deviant stimulus is 0.1
lambda.S=0.5*ones(1,300*length(stream));
lambda.D=0.5*ones(1,300*length(stream));
for i=1:length(stream)
   if(stream(i)=='S')
       lambda.S(300*(i-1)+1:300*(i-1)+50)=10;
       lambda.D(300*(i-1)+1:300*(i-1)+50)=2.5;
   elseif(stream(i)=='D')
       lambda.S(300*(i-1)+1:300*(i-1)+50)=2.5;
       lambda.D(300*(i-1)+1:300*(i-1)+50)=10;        
   end
end

%Stimulus generation
spike.ThS=binornd(ones(1,length(lambda.S)),lambda.S/1000);
spike.ThD=binornd(ones(1,length(lambda.D)),lambda.D/1000);

%Running the model with plasticity
[V,spike,w] = verP(spike,tau,w);

figure(9)
hold on
plot((1/60000):(1/60000):60,w.ThS.L4)
plot((1/60000):(1/60000):60,w.ThD.L4)
plot((1/60000):(1/60000):60,w.SP.L4)
legend('w.ThS->L4','w.ThD->L4','w.SP->L4')
xlabel('time(min)')
xlim([0 60])
ylabel('Synaptic weights')
title('4.Long term plasticity (p(S)=0.9 p(D)=0.1)')

%% 4.b.Running verNOP at points of interest
w1=w;                       %Make a copy of the weights
t=1e-3:1e-3:4.5;
temp.ThS.L4=[0.02 0.0326 0.0340 0.0349 0.0349 0.303];
temp.ThD.L4=[0.02 0.0741 0.1260 0.3011 0.4 0.4];
temp.SP.L4=[0.11 0.1051 0.1084 0.0352 0.0275 0.0006];
for i=1:6
    w1.ThS.L4=temp.ThS.L4(i);
    w1.ThD.L4=temp.ThD.L4(i);
    w1.SP.L4=temp.SP.L4(i);
    [~,spike1]=verNOP(spike1,tau,w1);       % Give it the earlier spike train used in verNOP
    figure(9+i)
    subplot(4,1,1)
    plot(t,spike1.ThS)
    ylabel('Spike.ThS')
    subplot(4,1,2)
    plot(t,spike1.ThD)
    ylabel('Spike.ThD')
    subplot(4,1,3)
    plot(t,spike1.SP)
    ylabel('Spike.SP')
    subplot(4,1,4)
    plot(t,spike1.L4)
    ylabel('Spike.L4')
    xlabel('t(s)')
    for j=1:4
       subplot(4,1,j)
       ylim([0 1])
    end
end
figure(10)
subplot(4,1,1)
title('4.Spiking activity - Original VerNOP')
figure(11)
subplot(4,1,1)
title('4.Spiking activity - Onset of D-weight increase (9.63 mins)')
figure(12)
subplot(4,1,1)
title('4.Spiking activity - SP-weight starts to fall (12.15 mins)')
figure(13)
subplot(4,1,1)
title('4.Spiking activity - SP-weight crosses S-weight (14.76 mins)')
figure(14)
subplot(4,1,1)
title('4.Spiking activity - Saturation of D-weight (15.23 mins)')
figure(15)
subplot(4,1,1)
title('4.Spiking activity - SP-L4 connection breaks (24.18 mins)')

%% 5.LTPlasticity with Equally likely S & D
stream=char('S'*ones(1,12000));
stream(binornd(ones(1,12000),0.5*ones(1,12000))==1)='D'; %Equal Probability of Standard and Deviant
lambda.S=0.5*ones(1,300*length(stream));
lambda.D=0.5*ones(1,300*length(stream));
for i=1:length(stream)
   if(stream(i)=='S')
       lambda.S(300*(i-1)+1:300*(i-1)+50)=10;
       lambda.D(300*(i-1)+1:300*(i-1)+50)=2.5;
   elseif(stream(i)=='D')
       lambda.S(300*(i-1)+1:300*(i-1)+50)=2.5;
       lambda.D(300*(i-1)+1:300*(i-1)+50)=10;        
   end
end

%Stimulus generation
spike.ThS=binornd(ones(1,length(lambda.S)),lambda.S/1000);
spike.ThD=binornd(ones(1,length(lambda.D)),lambda.D/1000);

%Running the model with plasticity
[V,spike,w] = verP(spike,tau,w);

figure(16)
hold on
plot((1/60000):(1/60000):60,w.ThS.L4)
plot((1/60000):(1/60000):60,w.ThD.L4)
plot((1/60000):(1/60000):60,w.SP.L4)
legend('w.ThS->L4','w.ThD->L4','w.SP->L4')
xlabel('time(min)')
xlim([0 60])
ylabel('Synaptic weights')
title('5.Long term plasticity (p(S)=p(D)=0.5)')
end

function [V,spike] = verNOP(spike,tau,w)
spike.SP=zeros(1,length(spike.ThS));
spike.L4=zeros(1,length(spike.ThS));
V.sp=zeros(1,length(spike.ThS));
V.l4=zeros(1,length(spike.ThS));
g.s=zeros(1,length(spike.ThS));
g.d=zeros(1,length(spike.ThS));
g.sp=zeros(1,length(spike.ThS));
x.es=zeros(1,length(spike.ThS));
x.rs=zeros(1,length(spike.ThS));
x.is=zeros(1,length(spike.ThS));
x.ed=zeros(1,length(spike.ThS));
x.rd=zeros(1,length(spike.ThS));
x.id=zeros(1,length(spike.ThS));
x.esp=zeros(1,length(spike.ThS));
x.rsp=zeros(1,length(spike.ThS));
x.isp=zeros(1,length(spike.ThS));
x.rs(1)=1;
x.rd(1)=1;
x.rsp(1)=1;
beta=5;

for t=2:length(spike.ThS)                   %time stamp 1 is initial condition
    % Update EPSP
    if(spike.ThS(t-1))
        g.s(t:end)=g.s(t:end)+exp(-((t:length(spike.ThS))-t)/tau.syn);        %if S spikes, add EPSP to g.s.
    end
    if(spike.ThD(t-1))
        g.d(t:end)=g.d(t:end)+exp(-((t:length(spike.ThS))-t)/tau.syn);        %if D spikes, add EPSP to g.d.
    end
    if(spike.SP(t-1))
        g.sp(t:end)=g.sp(t:end)+exp(-((t:length(spike.ThS))-t)/tau.syn);        %if SP spikes, add EPSP to g.d.
    end
    
    % Update membrane potential and account for leakage
    V.sp(t)=V.sp(t)+g.s(t-1)*w.ThS.SP*x.es(t-1)+g.d(t-1)*w.ThD.SP*x.ed(t-1);
    V.sp(t)=V.sp(t)-0.1*abs(V.sp(t-1));
    V.l4(t)=V.l4(t)+g.s(t-1)*w.ThS.L4*x.es(t-1)+g.d(t-1)*w.ThD.L4*x.ed(t-1)+g.sp(t-1)*w.SP.L4*x.esp(t-1);
    V.l4(t)=V.l4(t)-0.1*abs(V.l4(t-1));
    
    % Spiking activity and refractory modeling
    if(V.sp(t)>0.05)
        spike.SP(t)=1;
        V.sp(t+1:t+20)=V.sp(t)-beta*exp(-((t+1:t+20)-t-1)/tau.ref);
    end
    if(V.l4(t)>0.05)
        spike.L4(t)=1;
        V.l4(t+1:t+20)=V.l4(t)-beta*exp(-((t+1:t+20)-t-1)/tau.ref);
    end
    
    % Update neurotransmitter fractions
    [x.es(t),x.rs(t),x.is(t)]=STPlast(x.es(t-1),x.rs(t-1),x.is(t-1),tau.re,tau.ei.Th,tau.ir,spike.ThS(t-1));
    [x.ed(t),x.rd(t),x.id(t)]=STPlast(x.ed(t-1),x.rd(t-1),x.id(t-1),tau.re,tau.ei.Th,tau.ir,spike.ThD(t-1));
    [x.esp(t),x.rsp(t),x.isp(t)]=STPlast(x.esp(t-1),x.rsp(t-1),x.isp(t-1),tau.re,tau.ei.SP,tau.ir,spike.SP(t-1));
end

end

function [V,spike,w] = verP(spike,tau,w)
spike.SP=zeros(1,length(spike.ThS));
spike.L4=zeros(1,length(spike.ThS));
V.sp=zeros(1,length(spike.ThS));
V.l4=zeros(1,length(spike.ThS));
g.s=zeros(1,length(spike.ThS));
g.d=zeros(1,length(spike.ThS));
g.sp=zeros(1,length(spike.ThS));
x.es=zeros(1,length(spike.ThS));
x.rs=zeros(1,length(spike.ThS));
x.is=zeros(1,length(spike.ThS));
x.ed=zeros(1,length(spike.ThS));
x.rd=zeros(1,length(spike.ThS));
x.id=zeros(1,length(spike.ThS));
x.esp=zeros(1,length(spike.ThS));
x.rsp=zeros(1,length(spike.ThS));
x.isp=zeros(1,length(spike.ThS));
x.rs(1)=1;
x.rd(1)=1;
x.rsp(1)=1;
w.ThS.L4=zeros(1,length(spike.ThS));
w.ThD.L4=zeros(1,length(spike.ThS));
w.SP.L4=zeros(1,length(spike.ThS));
w.ThS.L4(1)=0.02;
w.ThD.L4(1)=0.02;
w.SP.L4(1)=0.11;
last.ThS=0;
last.ThD=0;
last.SP=0;
last.L4=0;
tau.ltp=13;
tau.ltd=20;
a.ltp=0.015;
a.ltd=0.021;
beta=5;

for t=2:length(spike.ThS)                   %time stamp 1 is initial condition
    % Update EPSP
    if(spike.ThS(t-1))
        g.s(t:end)=g.s(t:end)+exp(-((t:length(spike.ThS))-t)/tau.syn);        %if S spikes, add EPSP to g.s.
    end
    if(spike.ThD(t-1))
        g.d(t:end)=g.d(t:end)+exp(-((t:length(spike.ThS))-t)/tau.syn);        %if D spikes, add EPSP to g.d.
    end
    if(spike.SP(t-1))
        g.sp(t:end)=g.sp(t:end)+exp(-((t:length(spike.ThS))-t)/tau.syn);        %if SP spikes, add EPSP to g.d.
    end
    
    % Update membrane potential and account for leakage
    V.sp(t)=V.sp(t)+g.s(t-1)*w.ThS.SP*x.es(t-1)+g.d(t-1)*w.ThD.SP*x.ed(t-1);
    V.sp(t)=V.sp(t)-0.1*abs(V.sp(t-1));
    V.l4(t)=V.l4(t)+g.s(t-1)*w.ThS.L4(t-1)*x.es(t-1)+g.d(t-1)*w.ThD.L4(t-1)*x.ed(t-1)+g.sp(t-1)*w.SP.L4(t-1)*x.esp(t-1);
    V.l4(t)=V.l4(t)-0.1*abs(V.l4(t-1));
    
    % Spiking activity and refractory modeling
    if(V.sp(t)>0.05)
        spike.SP(t)=1;
        V.sp(t+1:t+20)=V.sp(t)-beta*exp(-((t+1:t+20)-t-1)/tau.ref);
    end
    if(V.l4(t)>0.05)
        spike.L4(t)=1;
        V.l4(t+1:t+20)=V.l4(t)-beta*exp(-((t+1:t+20)-t-1)/tau.ref);
    end
    
    % Update neurotransmitter fractions
    [x.es(t),x.rs(t),x.is(t)]=STPlast(x.es(t-1),x.rs(t-1),x.is(t-1),tau.re,tau.ei.Th,tau.ir,spike.ThS(t-1));
    [x.ed(t),x.rd(t),x.id(t)]=STPlast(x.ed(t-1),x.rd(t-1),x.id(t-1),tau.re,tau.ei.Th,tau.ir,spike.ThD(t-1));
    [x.esp(t),x.rsp(t),x.isp(t)]=STPlast(x.esp(t-1),x.rsp(t-1),x.isp(t-1),tau.re,tau.ei.SP,tau.ir,spike.SP(t-1));
    
    % Weight update
    w.ThS.L4(t)=w.ThS.L4(t-1);
    w.ThD.L4(t)=w.ThD.L4(t-1);
    w.SP.L4(t)=w.SP.L4(t-1);
    if(spike.L4(t))
        last.L4=t;
        if(last.ThS)
            w.ThS.L4(t)=w.ThS.L4(t-1)+w.ThS.L4(t-1)*a.ltp*exp(-(last.L4-last.ThS)/tau.ltp);
        end
        if(last.ThD)
            w.ThD.L4(t)=w.ThD.L4(t-1)+w.ThD.L4(t-1)*a.ltp*exp(-(last.L4-last.ThD)/tau.ltp);
        end
        if(last.SP)
            w.SP.L4(t)=w.SP.L4(t-1)+w.SP.L4(t-1)*a.ltp*exp(-(last.L4-last.SP)/tau.ltp);
        end
    end
    if(spike.ThS(t))
        last.ThS=t;
        if(last.L4)
            w.ThS.L4(t)=w.ThS.L4(t)-w.ThS.L4(t-1)*a.ltd*exp((last.L4-last.ThS)/tau.ltd);
        end
    end
    if(spike.ThD(t))
        last.ThD=t;
        if(last.L4)
            w.ThD.L4(t)=w.ThD.L4(t)-w.ThD.L4(t-1)*a.ltd*exp((last.L4-last.ThD)/tau.ltd);
        end
    end
    if(spike.SP(t))
        last.SP=t;
        if(last.L4)
            w.SP.L4(t)=w.SP.L4(t)-w.SP.L4(t-1)*a.ltd*exp((last.L4-last.SP)/tau.ltd);
        end
    end
    
    % Weight rescaling
    w.ThS.L4(t)=min(w.ThS.L4(t),0.4);
    w.ThS.L4(t)=max(w.ThS.L4(t),0.0001);
    w.ThD.L4(t)=min(w.ThD.L4(t),0.4);
    w.ThD.L4(t)=max(w.ThD.L4(t),0.0001);
    w.SP.L4(t)=min(w.SP.L4(t),0.11);
    w.SP.L4(t)=max(w.SP.L4(t),0.0001); 
end
end

function [x_e,x_r,x_i] = STPlast(x_eo,x_ro,x_io,t_re,t_ei,t_ir,spike)
    x_r=x_ro-spike*x_ro/t_re+x_io/t_ir;
    x_e=x_eo+spike*x_ro/t_re-x_eo/t_ei;
    x_i=x_io+x_eo/t_ei-x_io/t_ir;
end
