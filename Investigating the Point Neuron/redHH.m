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
