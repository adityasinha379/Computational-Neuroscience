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
