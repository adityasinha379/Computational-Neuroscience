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