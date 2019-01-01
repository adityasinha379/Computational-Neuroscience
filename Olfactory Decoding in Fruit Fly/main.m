clearvars; close all; clc;
rng(20184020);
%%
signal_generator = @(t, sample, sample_time,omega) ...
                    sum( diag(sample) * sinc( (repmat(t, size(sample_time))- ...
                    repmat(sample_time, size(t)))/pi*omega)/pi*omega);
                
M = 50; % sensory neurons
N = 20;  % projection neurons

h_OSN = @(t, hmax_OSN, tau_OSN) max(hmax_OSN*exp(-t/tau_OSN), 0).*(t>=0);
h_PN = @(t, hmax_PN, tau_PN) max(hmax_PN*(t/tau_PN).* exp(1-50*t/tau_PN), 0).*(t>=0);

dt     = 1e-5;                                      % set time step
t      = 0:dt:0.5;                                 % set simulation time
omega  = 2*pi*25;                                   % set cutoff frequency
T      = pi/omega;                                  % set sampling period
st     = (floor(t(1)/T)+1:floor(t(end)/T)+1)*T;     % set sample times
s      = rand(numel(st), 1);                    % set samples of the stimulus
u      = signal_generator(t, s, st', omega);        % generate the stimulus
u      = u/max(abs(u));                             % normalize the stimulus

%% OSN encoding
delta_l1 = 0.012*ones(M,1);
kappa_l1 = ones(M,1);
b_l1 = ones(M,1);

h_l1 = zeros(M,numel(t));

for i=1:M
    h_l1(i,:) = h_OSN(t,1+4*rand(1),(15+10*rand(1))*1e-3);
end

tk_OSN = OSN(t,u,h_l1,b_l1,kappa_l1,delta_l1);
spkavg_l1 = sum(cellfun(@numel,tk_OSN))/M;

spike_l1 = zeros(M,length(t));
for i=1:M
    spike_l1(i,round(tk_OSN{i}/dt+1))=1;
end

%% Relevant plots
figure(1)
subplot(311)
plot(t,u,'r')
title('Input signal');xlabel('t(s)');
subplot(312)
plot(t,h_l1(1,:),'r')
title('Filter h_1');xlabel('t(s)');
subplot(313)
temp = conv(u,h_l1(1,:))*dt;
plot(t,temp(1:length(t)),'r')
title('Filtered signal');xlabel('t(s)')

figure(2)
utils.plotRaster(t,spike_l1)
title('OSN Spiking Raster plot')
xlabel('t(s)')

figure(3)
hold on
stem(1:M,cellfun(@numel,tk_OSN)','r');
line([0 M],[50 50],'color','b');
xlabel('OSN number')
ylabel('Spike Rate (Hz)')
title('OSN Spiking')
legend('Spike rates','Nyquist rate')

%% Glomerulus
h_l2 = cell(N,M);

for i=1:N
    for j=1:M
        h_l2{i,j} = h_PN(t,1+4*rand(1),(350+100*rand(1))*1e-3);
    end
end

v_PN = zeros(N,length(t));

for i=1:N
    for j=1:M
        temp = conv(spike_l1(j,:),h_l2{i,j});
        v_PN(i,:) = v_PN(i,:) + temp(1:length(t));
    end
end

%% Projection Neurons
delta_l2 = 0.04*ones(N,1);
kappa_l2 = ones(N,1);
b_l2 = 1.7*ones(N,1);

tk_PN = PN(t,v_PN,b_l2,kappa_l2,delta_l2);
spkavg_l2 = sum(cellfun(@numel,tk_PN))/N;

spike_l2 = zeros(N,length(t));
for i=1:N
    spike_l2(i,round(tk_PN{i}/dt+1))=1;
end
%% Relevant Plots
figure(4)
hold on
for i=1:N
    for j=1:M
        plot(t,h_l2{i,j})
    end
end
title('Glomerular filters h_n_m');
xlabel('t(s)')

figure(5)
hold on
for i=1:N
   plot(t,v_PN(i,:))
end
title('Glomerular outputs v_1_n');
xlabel('t(s)')

figure(6)
utils.plotRaster(t,spike_l2)
title('OSN Spiking Raster plot')
xlabel('t(s)')

figure(7)
hold on
stem(1:N,cellfun(@numel,tk_PN)','r');
line([0 N],[100 100],'color','b');
xlabel('PN number')
ylabel('Spike Rate (Hz)')
title('PN Spiking')
legend('Spike rates','Nyquist rate')

%% OSN Decoding & Recovery
dim = 0;
for i=1:M
   dim = dim+length(tk_OSN{i})-1;
end

temp = zeros(1,M+1);
for i=2:M+1         %Find out the cumulative lengths of tk
    temp(i)=temp(i-1)+length(tk_OSN{i-1})-1;
end

q_l1 = zeros(dim,1);
G_l1 = zeros(dim);

for i=1:M
    q_l1(temp(i)+1:temp(i)+length(tk_OSN{i})-1) = kappa_l1(i)*delta_l1(i) - b_l1(i)*diff(tk_OSN{i})';
end
%%
t_long = -0.5:dt:0.5;

for j=1:M
    cum_hg = dt*cumtrapz(dt*conv(h_l1(j,:),sinc(t_long/T)/T));
    sk = (tk_OSN{j}(1:end-1)+tk_OSN{j}(2:end))/2;
    for i =1:M
        for k = 1:length(tk_OSN{j})-1
            G_l1(temp(i)+1:temp(i)+length(tk_OSN{i})-1,temp(j)+k) = diff(cum_hg(round((tk_OSN{i} -sk(k) - t_long(1))/dt)+1));
        end
    end
end
%%
% Reconstruction
u_rec = zeros(size(t));

c = pinv(G_l1,1e-3)*q_l1;
for j=1:M
    sk = (tk_OSN{j}(1:end-1)+tk_OSN{j}(2:end))/2;
    u_rec = u_rec + signal_generator(t,c(temp(j)+1:temp(j+1)),sk',omega);
end
%% Relevant plots
figure('Name','Decoding OSN Stage','Color','w','Position',[10 30 1200 400]);
subplot(211);plot(1e3*t,u,'r',1e3*t,u_rec,'--b','LineWidth',2);
title('OSN recovery (u(t))'); xlabel('t(ms)')
legend('Input stimulus','Recovered stimulus','Location','SouthEast')
set(subplot(212),'ylim',[-120 50],'title',text('string','Recovery Error (dB)'),...
        'xlabel',text('string','t(ms)'))
    hold on;plot(1e3*t,mag2db(sqrt((u-u_rec).^2)),'-k','LineWidth',2);

%% PN Decoding & Recovery
v_PN_rec = zeros(N,numel(t));
omega_new = 2*pi*100;
T_new = pi/omega_new;

for i=1:N
    q_l2 = kappa_l2(i)*delta_l2(i) - b_l2(i)*diff(tk_PN{i}');
    sk = (tk_PN{i}(2:end) + tk_PN{i}(1:end-1))/2;
    G_l2 = zeros(numel(tk_PN{i})-1);
    for j = 1:numel(sk)
        temp = dt*cumtrapz( sinc((t-sk(j))/T_new)/T_new );
        G_l2(:,j) = diff( temp( round( (tk_PN{i}-t(1))/dt )+1 ) );
    end
    
    c = pinv(G_l2,1e-6)*q_l2;

    v_PN_rec(i,:) = signal_generator(t,c,sk',omega_new);
end

%%
figure('Name','Decoding PN Stage','Color','w','Position',[10 30 1200 400]);
    set(subplot(2,1,1),'title',text('string',{'PN stage decoding of v_1_1',}));
    hold on;plot(1e3*t,v_PN(1,:),'r',1e3*t,v_PN_rec(1,:),'--b','LineWidth',2); 
    legend('Stimulus','Recovered Stimulus','Location','SouthEast');    
    set(subplot(2,1,2),'title',text('string','Recovery Error (dB)'),...
        'xlabel',text('string','t(ms)'))
    hold on;plot(1e3*t,mag2db(sqrt((v_PN(1,:)-v_PN_rec(1,:)).^2)),'-k','LineWidth',2);
    
%% Parallel Filters
u_fft = fft(u,numel(t))/numel(t);
v_PN_fft = zeros(N,numel(t));
h_rec = zeros(N,numel(t));
f = 0:1/dt/numel(t):1/dt/2-1/dt/numel(t);

for i=1:1
   v_PN_fft(i,:) = fft(v_PN(i,:),numel(t))/numel(t);
   h_fft = v_PN_fft(i,:)./u_fft;
   h_rec(i,:) = ifft(h_fft);
end

figure;
subplot(311)
plot(f(1:100),abs(v_PN_fft(1,1:length(f(1:100)))),'r');
xlabel('f (Hz)'); title('FFT of v_1_1');
subplot(312)
plot(f(1:100),abs(u_fft(1:length(f(1:100)))),'r');
xlabel('f (Hz)'); title('FFT of u');
subplot(313)
plot(f(1:100),abs(h_fft(1:length(f(1:100)))),'r');
xlabel('f (Hz)'); title('FFT of h_r_e_c');

%% End-to-End Recovery
dim = 0;
for i=1:N
   dim = dim+length(tk_PN{i})-1;
end

temp = zeros(1,N+1);
for i=2:N+1         %Find out the cumulative lengths of tk
    temp(i)=temp(i-1)+length(tk_PN{i-1})-1;
end

q = zeros(dim,1);
G = zeros(dim);

for i=1:N
    q(temp(i)+1:temp(i)+length(tk_PN{i})-1) = kappa_l2(i)*delta_l2(i) - b_l2(i)*diff(tk_PN{i})';
end
%%
t_long = -0.5:dt:0.5;

for j=1:N
    cum_hg = dt*cumtrapz(dt*conv(h_rec(j,:),sinc(t_long/T)/T));
    sk = (tk_PN{j}(1:end-1)+tk_PN{j}(2:end))/2;
    for i =1:N
        for k = 1:length(tk_OSN{j})-1
            G(temp(i)+1:temp(i)+length(tk_PN{i})-1,temp(j)+k) = diff(cum_hg(round((tk_PN{i} -sk(k) - t_long(1))/dt)+1));
        end
    end
end
%%
% Reconstruction
u_rec_e2e = zeros(size(t));

c = pinv(G,1e-3)*q;
for j=1:N
    sk = (tk_PN{j}(1:end-1)+tk_PN{j}(2:end))/2;
    u_rec_e2e = u_rec_e2e + signal_generator(t,c(temp(j)+1:temp(j+1)),sk',omega)/2e2;
end

% Relevant plots
figure('Name','End-to-End Recovery','Color','w','Position',[10 30 1200 400]);
subplot(311);plot(1e3*t,u,'r','LineWidth',2);
title('OSN input stimulus (u(t))'); xlabel('t(ms)')
subplot(312);plot(1e3*t,u_rec_e2e,'--b','LineWidth',2);
title('End-to-End Recovered stimulus'); xlabel('t(ms)')
set(subplot(313),'ylim',[-120 50],'title',text('string','Recovery Error (dB)'),...
        'xlabel',text('string','t(ms)'))
    hold on;plot(1e3*t,mag2db(sqrt((u-u_rec_e2e).^2)),'-k','LineWidth',2);