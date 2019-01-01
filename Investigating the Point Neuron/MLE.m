function dxdt = MLE(t,x)      %MLE equation description solved by ODE15s
global Gca Gk Gl Eca Ek El phi V1 V2 V3 V4 C Iext;

V=x(1);
w=x(2);

m_inf=0.5*(1+tanh((V-V1)/V2));
w_inf=0.5*(1+tanh((V-V3)/V4));
tau_w=1/cosh((V-V3)/(2*V4));

dxdt=zeros(2,1);
dxdt(1)=(1/C)*(Iext-Gca*m_inf*(V-Eca)-Gk*w*(V-Ek)-Gl*(V-El));
dxdt(2)=phi*(w_inf-w)/tau_w;
end
