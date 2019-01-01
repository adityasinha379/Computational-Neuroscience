function HHpart
close all;
clear;
%
%% 12. Hodgkin Huxley Model
global Gk Gna Gl Ek Ena El C Iext Veq meq heq neq f;
C = 1; Iext = 0;
Gna=120; Gk=36; Gl=0.3; Ena=55; Ek=-72; El=-49.4038;
%}
%
%% 13. Determining El for Vr = -60mV
syms V m h n;
V_nc=(Iext-Gk*(n^4)*(V-Ek)-Gna*(m^3)*h*(V-Ena)-Gl*(V-El))==0;
m_nc=(((-0.1) * (35+V) / (exp(-0.1*(35+V)) - 1))*(1-m)-(4 * exp(-(60+V)/18))*(m))==0;
h_nc=((0.07 * exp(-(60+V)/20))*(1-h)-(1 / (exp(-(30+V)/10) + 1))*(h))==0;
n_nc=(((-0.01) * (50+V) / (exp(-0.1*(50+V)) - 1))*(1-n)-(0.125 * exp(-(60+V)/80))*(n))==0;
warning('off','symbolic:solve:FallbackToNumerical')
[Veq,meq,heq,neq]=vpasolve([V_nc,m_nc,h_nc,n_nc],[V,m,h,n]);    
Veq = double(Veq); meq = double(meq); heq = double(heq); neq = double(neq);
Elq = (-60)+((Gk*(neq^4)*(-60-Ek)+Gna*(meq^3)*heq*(-60-Ena))/Gl);
%}
%
%% 14. Determining Action Potential Threshold  
figure(15)
 tspan=[0 100];
Vs=[-53.495:0.001:-53.49];
for i=1:length(Vs)
x0=[Vs(i);meq;heq;neq];
[t,x]=ode15s(@HH,tspan,x0);
plot(t,x(:,1));
hold on
end
h=get(gca,'Children');
legend([h(6) h(5) h(4) h(3) h(2) h(1)],['V=' num2str(Vs(1))],['V=' num2str(Vs(2))],...
    ['V=' num2str(Vs(3))],['V=' num2str(Vs(4))],['V=' num2str(Vs(5))],['V=' num2str(Vs(6))])
title('Determining Action Potential Threshold')
Vthr = -53.493;
%}
%
%% 15. Effect of External current
figure(16)
I=[8:12];
for i=1:length(I)
syms V m h n;
V_nc=(I(i)-Gk*(n^4)*(V-Ek)-Gna*(m^3)*h*(V-Ena)-Gl*(V-El))==0;
m_nc=(((-0.1) * (35+V) / (exp(-0.1*(35+V)) - 1))*(1-m)-(4 * exp(-(60+V)/18))*(m))==0;
h_nc=((0.07 * exp(-(60+V)/20))*(1-h)-(1 / (exp(-(30+V)/10) + 1))*(h))==0;
n_nc=(((-0.01) * (50+V) / (exp(-0.1*(50+V)) - 1))*(1-n)-(0.125 * exp(-(60+V)/80))*(n))==0;
warning('off','symbolic:solve:FallbackToNumerical')
[Veq,meq,heq,neq]=vpasolve([V_nc,m_nc,h_nc,n_nc],[V,m,h,n]);    
Veq = double(Veq);
meq = double(meq);
heq = double(heq);
neq = double(neq);   

jac=jacobian([(1/C)*(Iext-Gk*(n^4)*(V-Ek)-Gna*(m^3)*h*(V-Ena)-Gl*(V-El)),...
    ((-0.1) * (35+V) / (exp(-0.1*(35+V)) - 1))*(1-m)-(4 * exp(-(60+V)/18))*(m)...
    (0.07 * exp(-(60+V)/20))*(1-h)-(1 / (exp(-(30+V)/10) + 1))*(h)...
    ((-0.01) * (50+V) / (exp(-0.1*(50+V)) - 1))*(1-n)-(0.125 * exp(-(60+V)/80))*(n)], [V,m,h,n]);
jac=double(subs(jac,[V,m,h,n],[Veq,meq,heq,neq]));
disp([' Iext = ' num2str(I(i)) '   Veq = ' num2str(Veq) '   meq = ' num2str(meq)...
        '   heq = ' num2str(heq) '   neq = ' num2str(neq)]);
eigenvalues=eig(jac)
x0=[Veq;meq+0.01;heq;neq];
[t,x]=ode15s(@HH,tspan,x0);
plot(t,x(:,1));
hold on
end
h=get(gca,'Children');
legend([h(5) h(4) h(3) h(2) h(1)],['I_e_x_t= ' num2str(I(1))],['I_e_x_t= ' num2str(I(2))],...
    ['I_e_x_t= ' num2str(I(3))],['I_e_x_t= ' num2str(I(4))],['I_e_x_t= ' num2str(I(5))])
title('Action Potentials for 8\muA/cm^2 < I_e_x_t < 12\muA/cm^2')
%}
%
%% 16. Action Potentials in case of myotonia
figure(17)
% fi=[0 0.1 0.17 0.2];
% for i=1:length(fi)
syms V m h n;
% f=fi(i);
f=0.1;
V_nc=(Iext-Gk*(n^4)*(V-Ek)-Gna*(1-f)*(m^3)*h*(V-Ena)-Gna*f*(m^3)*(V-Ena)-Gl*(V-El))==0;
m_nc=(((-0.1) * (35+V) / (exp(-0.1*(35+V)) - 1))*(1-m)-(4 * exp(-(60+V)/18))*(m))==0;
h_nc=((0.07 * exp(-(60+V)/20))*(1-h)-(1 / (exp(-(30+V)/10) + 1))*(h))==0;
n_nc=(((-0.01) * (50+V) / (exp(-0.1*(50+V)) - 1))*(1-n)-(0.125 * exp(-(60+V)/80))*(n))==0;
warning('off','symbolic:solve:FallbackToNumerical')
[Veq,meq,heq,neq]=vpasolve([V_nc,m_nc,h_nc,n_nc],[V,m,h,n]);    
Veq = double(Veq);
meq = double(meq);
heq = double(heq);
neq = double(neq);   

Vs=[-65:1:-60];
for i=1:length(Vs)
x0=[Vs(i);meq;heq;neq];
[t,x]=ode15s(@HHmyotonia,tspan,x0);
plot(t,x(:,1));
hold on
end
h=get(gca,'Children');
legend([h(6) h(5) h(4) h(3) h(2) h(1)],['V=' num2str(Vs(1))],['V=' num2str(Vs(2))],...
    ['V=' num2str(Vs(3))],['V=' num2str(Vs(4))],['V=' num2str(Vs(5))],['V=' num2str(Vs(6))])
title('Action Potential')
%}
%
%% 17. Action Potential with Reduced HH system 
figure(18)
tspan=[0 100];
Vs=[-55:1:-50];
for i=1:length(Vs)
x0=[Vs(i);neq];
[t,x]=ode15s(@redHH,tspan,x0);
plot(t,x(:,1));
hold on
end
h=get(gca,'Children');
legend([h(6) h(5) h(4) h(3) h(2) h(1)],['V=' num2str(Vs(1))],['V=' num2str(Vs(2))],...
    ['V=' num2str(Vs(3))],['V=' num2str(Vs(4))],['V=' num2str(Vs(5))],['V=' num2str(Vs(6))])
title('Action Potential using reduced HH')

figure(19)
tspan=[0 100];
Vs=[-55:1:-50];
for i=1:length(Vs)
x0=[Vs(i);meq;heq;neq];
[t,x]=ode15s(@HH,tspan,x0);
plot(t,x(:,1));
hold on
end
h=get(gca,'Children');
legend([h(6) h(5) h(4) h(3) h(2) h(1)],['V=' num2str(Vs(1))],['V=' num2str(Vs(2))],...
    ['V=' num2str(Vs(3))],['V=' num2str(Vs(4))],['V=' num2str(Vs(5))],['V=' num2str(Vs(6))])
title('Action Potential using HH')
%}
%
%% 18. Study of Bifurcation
fi = [0.02:0.005:0.10 0.11:0.1:0.4];
for i=1:length(fi)
figure(20)
syms V
f=fi(i);
m_inf =@(V) ((-0.1) * (35+V) / (exp(-0.1*(35+V)) - 1))/(((-0.1) * (35+V) / (exp(-0.1*(35+V)) - 1))+(4 * exp(-(60+V)/18)));
h_inf =@(V) (0.07 * exp(-(60+V)/20))/((0.07 * exp(-(60+V)/20))+(1 / (exp(-(30+V)/10) + 1)));
Vnc =@(V) ((Iext-Gna*(1-f)*(m_inf^3)*h_inf*(V-Ena)-Gna*(f)*(m_inf^3)*(V-Ena)...
    -Gl*(V-El))/(Gk*(V-Ek))^0.25);
n_nc=@(V) (((-0.01)*(50+V)/(exp(-0.1*(50+V))-1))/(((-0.01)*(50+V)/(exp(-0.1*(50+V))-1))+0.125*exp(-(60+V)/80)));

fplot(@(V) (V_nc(V)^4)*100,[-84 120 0 100]);
hold on
fplot(@(V) n_nc(V)*100,[-84 120 0 100]);
hold on
end
xlabel('V(mV)')
ylabel('w*100')
title('Phase-Plane Plot')
%}
end

% 12) HH equations
function dxdt = HH(t,x)              
global Gk Gna Gl Ek Ena El C Iext am ah an bm bh bn ;
V=x(1);
m=x(2);
h=x(3);
n=x(4);
am = (-0.1) * (35+V) / (exp(-0.1*(35+V)) - 1);
bm = 4 * exp(-(60+V)/18);
ah = 0.07 * exp(-(60+V)/20);
bh = 1 / (exp(-(30+V)/10) + 1);
an = (-0.01) * (50+V) / (exp(-0.1*(50+V)) - 1);
bn = 0.125 * exp(-(60+V)/80);

dxdt=zeros(4,1);
dxdt(1)=(1/C)*(Iext-Gk*(n^4)*(V-Ek)-Gna*(m^3)*h*(V-Ena)-Gl*(V-El));
dxdt(2)=am*(1-m)-bm*(m);
dxdt(3)=ah*(1-h)-bh*(h);
dxdt(4)=an*(1-n)-bn*(n);
end

% 16. HH myotonia equations
function dydt = HHmyotonia(t,x)      
global Gk Gna Gl Ek Ena El C Iext f;
V=x(1);
m=x(2);
h=x(3);
n=x(4);
am = (-0.1) * (35+V) / (exp(-0.1*(35+V)) - 1);
bm = 4 * exp(-(60+V)/18);
ah = 0.07 * exp(-(60+V)/20);
bh = 1 / (exp(-(30+V)/10) + 1);
an = (-0.01) * (50+V) / (exp(-0.1*(50+V)) - 1);
bn = 0.125 * exp(-(60+V)/80);

dydt=zeros(4,1);
dydt(1)=(1/C)*(Iext-Gk*(n^4)*(V-Ek)-Gna*(1-f)*(m^3)*h*(V-Ena)-Gna*(f)*(m^3)*(V-Ena)-Gl*(V-El));
dydt(2)=am*(1-m)-bm*(m);
dydt(3)=ah*(1-h)-bh*(h);
dydt(4)=an*(1-n)-bn*(n);
end

% 17. Reduced HH equations
function dxdt = redHH(t,x)           
global Gk Gna Gl Ek Ena El C Iext ;
V=x(1);
n=x(2);
am = (-0.1) * (35+V) / (exp(-0.1*(35+V)) - 1);
bm = 4 * exp(-(60+V)/18);
ah = 0.07 * exp(-(60+V)/20);
bh = 1 / (exp(-(30+V)/10) + 1);
an = (-0.01) * (50+V) / (exp(-0.1*(50+V)) - 1);
bn = 0.125 * exp(-(60+V)/80);
m_inf = am/(am+bm);
h_inf = ah/(ah+bh);

dxdt=zeros(2,1);
dxdt(1)=(1/C)*(Iext-Gk*(n^4)*(V-Ek)-Gna*(m_inf^3)*h_inf*(V-Ena)-Gl*(V-El));
dxdt(2)=an*(1-n)-bn*(n);
end
