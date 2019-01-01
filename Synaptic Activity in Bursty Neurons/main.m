clearvars; close all; clc;
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
rng(20184020);
fontsize = 12;
linewidth = 1.5;
%% 1. Simulating the Bursty model
spike_detect_online = @(V,threshold) (V(2)>V(1)) & (V(2) > V(3)) & (V(2) > threshold);
spike_detect = @(V,threshold) [false,(V(2:end-1)>V(1:end-2)) & (V(2:end-1) > V(3:end)) & (V(2:end-1) > threshold),false];

dt = 1e-6; % Set time step
t = -0.05:dt:0.4; % Create Time Signal
Iext = 35.6*ones(size(t)); % Create fixed input current

global syn;
syn.dualsyn=0;                  % Controls whether we call one or two bursty neurons
out=bursty(t,Iext,[-36;0;1]);

figure(1)
view(axes(),3)
title('Dynamics of Bursty Neuron');
xlim([0 1]); ylim([0 1]); zlim([-50 10]);
xlabel('m'); ylabel('c'); zlabel('V (mV)');
hold on
plot3(out(2,:),out(3,:),out(1,:),'r');
temp=[];
for i=1:200:length(t)
    delete(temp);
    temp=plot3(out(2,i),out(3,i),out(1,i),'ko');
    set(temp,'markerfacecolor','g');
    pause(0.01);
end

figure(2)
suptitle('Bursty neuron time traces')
for i=1:size(out,1)
   subplot(size(out,1),1,i)
   xlabel('t(s)');
   plot(t,out(i,:))
end
subplot(size(out,1),1,1);
ylabel('V (mV)');
subplot(size(out,1),1,2);
ylabel('m');
subplot(size(out,1),1,3);
ylabel('c');
V.bursty=out(1,:); x.bursty=out(2:3,:); out=[];


[V.hh,x.hh]=hodgkin_huxley(t,7*ones(size(t)));

figure(3)
view(axes(),3)
title('Dynamics of HH Neuron');
xlim([0 1]); ylim([0 1]); zlim([-20 120]);
xlabel('n'); ylabel('h'); zlabel('V (mV)');
temp=[];
hold on
plot3(x.hh(2,:),x.hh(3,:),V.hh,'r');
for i=1:200:length(t)
    delete(temp);
    temp=plot3(x.hh(2,i),x.hh(3,i),V.hh(i),'ko');
    set(temp,'markerfacecolor','g');
    pause(0.01);
end

%% 2. Connecting Bursty model to Rinzel model
syn.E_pre=2; syn.E_rev=-70*ones(1,10); syn.k_pre=0.22; syn.alpha=5000; syn.beta=0.18;
syn.Gsyn=0.03; syn.Tmax=0.002;

%% a. Inhibitory synapses for different Gsyn values
Iext=15*ones(size(t));
Gsyn=[0.03 .15 .5 1 2];

[V.rinz,~]=hh_rinzel(t,Iext,[]);

figure(4)
suptitle('a. Rinzel responses for different max conductances (inhibitory)');
for i=1:length(Gsyn)
    subplot(length(Gsyn),1,i);
    hold on
    plot(t,V.rinz);
    syn.Gsyn=Gsyn(i);
    [V.rinz_syn,~]=hh_rinzel(t,Iext,repmat(V.bursty,size(syn.E_rev,1)));
    plot(t,V.rinz_syn)
    xlabel('t(s)')
    ylabel('V(mV)')
    legend('Independent Rinzel neuron',['G_s_y_n = ' num2str(Gsyn(i))])
end

%% b. Excitatory synapses for different Gsyn values
syn.E_rev=0*ones(1,10);

figure(5)
suptitle('b. Rinzel responses for different max conductances (excitatory)');
for i=1:length(Gsyn)
    subplot(length(Gsyn),1,i);
    hold on
    plot(t,V.rinz);
    syn.Gsyn=Gsyn(i);
    [V.rinz_syn,~]=hh_rinzel(t,Iext,repmat(V.bursty,size(syn.E_rev,1)));
    plot(t,V.rinz_syn)
    xlabel('t(s)')
    ylabel('V(mV)')
    legend('Independent Rinzel neuron',['G_s_y_n = ' num2str(Gsyn(i))]);
end

%% c. Mixture of excitatory and inhibitory synapses
syn.Gsyn=0.03;

figure(6)
suptitle('c. Rinzel responses for mixture of exctitatory and inhibitory inputs');

for i=1:5
    subplot(5,1,i);
    hold on
    syn.E_rev=-70*ones(1,10).*randi([0 1],1,10);
    temp=sum(-syn.E_rev/70);
    plot(t,V.rinz);
    [V.rinz_syn,~]=hh_rinzel(t,Iext,repmat(V.bursty,size(syn.E_rev,1)));
    plot(t,V.rinz_syn)
    xlabel('t(s)')
    ylabel('V(mV)')
    title([num2str(10-temp) ' excitatory, ' num2str(temp) ' inhibitory']);
    legend('Independent Rinzel neuron',['G_s_y_n = ' num2str(syn.Gsyn)]);
end

%% d. Excitatory synapses with different alpha values
syn.E_rev=0*ones(1,10);
alpha=[1000 3000 5000 7000 9000];

figure(7)
suptitle('d. Rinzel responses - excitatory for different \alpha');

for i=1:length(alpha)
    subplot(length(alpha),1,i);
    hold on
    plot(t,V.rinz);
    syn.alpha=alpha(i);
    [V.rinz_syn,~]=hh_rinzel(t,Iext,repmat(V.bursty,size(syn.E_rev,1)));
    plot(t,V.rinz_syn)
    xlabel('t(s)')
    ylabel('V(mV)')
    legend('Independent Rinzel neuron',['$\alpha$ = ' num2str(syn.alpha)]);
end

%% e. Time varying alpha
syn.alpha=5^3*ones(size(t));
[~,temp] = min(abs(t));
syn.alpha(temp+1:temp+round(0.2/dt))=8^4;

figure(8)
suptitle('e. Rinzel responses - excitatory for time-varying $\alpha$');
subplot(211)
plot(t,syn.alpha);
subplot(212)
hold on
plot(t,V.rinz);
[V.rinz_syn,~]=hh_rinzel(t,Iext,repmat(V.bursty,size(syn.E_rev,1)));
plot(t,V.rinz_syn)
xlabel('t(s)')
ylabel('V(mV)')
legend('Independent Rinzel neuron','Rinzel w Bursty input');

%% 3. Connecting two bursty neurons together
% a. Maximum conductance of two inhibitory synapses
syn.aplha=5000;
Iext=35.6*ones(size(t));
syn.dualsyn=1;
Gsyn=[0.01 0.03 0.1];

figure(9)
suptitle('Interconnected bursty neurons: Time traces');
subplot(3,1,1);ylabel('V (mV)');
hold on; plot(t,V.bursty);
subplot(3,1,2);ylabel('m');title('Neuron 1/2');
hold on; plot(t,x.bursty(1,:));
subplot(3,1,3);ylabel('c');
hold on; plot(t,x.bursty(2,:));

legendinf=cell(1,length(Gsyn)+1);
legendinf{1}='Independent bursty neuron';
for j=1:length(Gsyn)
    syn.E_rev=[-70;-70];
    syn.Gsyn=Gsyn(j);
    out=bursty(t,Iext,[-36;0;1]);
    for i=1:3
       subplot(3,1,i)
       xlabel('t(s)');
       plot(t,out(i,:))
    end
    legendinf{j+1}=['G_s_y_n = ' num2str(Gsyn(j))];
end
legend(legendinf,'Position',[0.85 0.32 0.1 0.1]);

%% b. Combination of excitatory and inhibitory
syn.Gsyn=0.03;
syn.dualsyn=1;
temp=[0 0; 0 1; 1 1];

figure(10)
suptitle('Interconnected bursty neurons: Time traces');
subplot(3,2,1);title('Neuron 1');
subplot(3,2,2);title('Neuron 2');
for j=1:size(temp,1)*size(temp,2)
    subplot(3,2,j)
    hold on; plot(t,V.bursty);
    ylabel('V (mV)');
end

for j=1:size(temp,1)
    syn.E_rev=-70*temp(j,:);
    out=bursty(t,Iext,[-36;0;1]);
    for i=1:2
       subplot(3,2,2*j-2+i)
       xlabel('t(s)');
       plot(t,out(3*i-2,:))
    end
end
for j=1:2
   subplot(3,2,j)
   legend('Independent bursty neuron', 'Both excitatory (II)') 
end
for j=3:4
   subplot(3,2,j)
   legend('Independent bursty neuron', 'Mixed (EI)')
end
for j=5:6
   subplot(3,2,j)
   legend('Independent bursty neuron', 'Both inhibitory (EE))')
end

%% c. Rate of synapse alpha
alpha = [1000 5000 9000];
syn.dualsyn=1;

figure(11)
suptitle('Interconnected bursty neurons: Time traces');
subplot(3,1,1);ylabel('V (mV)');
hold on; plot(t,V.bursty);
subplot(3,1,2);ylabel('m');title('Neuron 1/2');
hold on; plot(t,x.bursty(1,:));
subplot(3,1,3);ylabel('c');
hold on; plot(t,x.bursty(2,:));

legendinf=cell(1,length(alpha)+1);
legendinf{1}='Independent bursty neuron';
for j=1:length(alpha)
    syn.E_rev=[-70;-70];
    syn.alpha=alpha(j);
    out=bursty(t,Iext,[-36;0;1]);
    for i=1:3
       subplot(3,1,i)
       xlabel('t(s)');
       plot(t,out(i,:))
    end
    legendinf{j+1}=['\alpha = ' num2str(alpha(j))];
end
legend(legendinf,'Position',[0.85 0.32 0.1 0.1]);

%% Creation of degenerate responses
syn.dualsyn=1;
syn.E_rev=[-70;-70];
% alpha and Gsyn
syn.alpha=5000;
syn.Gsyn=0.03;

out=bursty(t,Iext,[-36;0;1]);
figure(12)
title('Interconnected bursty neurons: Degenerate responses');
xlabel('t(s)');ylabel('V (mV)');
hold on;
plot(t,out(1,:));
legendinf=cell(1,2);
legendinf{1}=(['G_s_y_n = ' num2str(syn.Gsyn) ', \alpha = ' num2str(syn.alpha)]);

syn.E_rev=[-70;-70];
syn.alpha=1000;
syn.Gsyn=0.06;
out=bursty(t,Iext,[-36;0;1]);
plot(t,out(1,:));
legendinf{2}=(['G_s_y_n = ' num2str(syn.Gsyn) ', \alpha = ' num2str(syn.alpha)]);
legend(legendinf);