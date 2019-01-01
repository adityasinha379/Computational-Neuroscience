function x = bursty(t,Iext,init)

global syn;
dt=1e3*diff(t(1:2));
n1.C=1.2;
n1.del=0.052;
n1.eps=4.9;

n1.g=[4.4 9.0 .19 2.0];    % Ca K KS L
n1.E=[120 -80 -80 -60];
n1.Vth=[-1.2 2.0 -27];
n1.ko=[1/18 0.1 0.8];
n1.Isyn=0;
n1.s=0;
if(syn.dualsyn==1)
    n2=n1;
    E_rev=syn.E_rev;
end

x=zeros(3+3*syn.dualsyn,length(t));   % V m c
if(nargin<3)
   init=[0;0;1]; 
end
x(:,1)=repmat(init,[1+syn.dualsyn,1]);
dx=zeros(3+3*syn.dualsyn,1);

for i=2:length(t)
    if(syn.dualsyn==1)
        syn.E_rev=E_rev(1);
        [n1.Isyn,n1.s] = synapse(n1.s,x(4,i-1),x(1,i-1));   % n2->n1
        syn.E_rev=E_rev(2);
        [n2.Isyn,n2.s] = synapse(n2.s,x(1,i-1),x(4,i-1));   % n1->n2
    end
    dx(1:3) = bursty_ode(x(1:3,i-1), Iext(i)-n1.Isyn, n1);
    if(syn.dualsyn==1)
        dx(4:6) = bursty_ode(x(4:6,i-1), Iext(i)-n2.Isyn, n2);
    end
    x(:,i) = x(:,i-1) + dt*dx;
end

end
