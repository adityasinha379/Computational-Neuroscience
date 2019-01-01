function Project2
close all;

%
%% 2. Nullclines, equilibrium point, quiver plot
global Gca Gk Gl Eca Ek El phi V1 V2 V3 V4 C Iext;

Gca=4.4; Gk=8.0; Gl=2; Eca=120; Ek=-84; El=-60;             %values of part 2
phi=0.02; V1=-1.2; V2=18; V3=2; V4=30; C=20; Iext=0;

figure(1)
hold on

syms V w;
V_nc=(Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gk*w*(V-Ek)-Gl*(V-El))==0;
w_nc=(0.5*(1+tanh((V-V3)/V4))-w)==0;
warning('off','symbolic:solve:FallbackToNumerical')
[Veq,weq]=solve([V_nc,w_nc],[V,w]);    % 1st method for eqm pt - intersection of nullclines
Veq1=Veq;weq1=weq;
Veq=double(Veq);
weq=double(weq);
h=plot(Veq,weq*100,'+');
set(h,'linewidth',2)
text(Veq,weq*100,['  \leftarrow(' num2str(round(Veq,5)) ',' num2str(round(weq,5)) ')'])

V_nc=@(V) (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gl*(V-El))/(Gk*(V-Ek));  % for nullclines
w_nc=@(V) (0.5*(1+tanh((V-V3)/V4)));
fplot(@(V) V_nc(V)*100,[-84 120 0 100]);
fplot(@(V) w_nc(V)*100,[-84 120 0 100]);
xlabel('V(mV)')
ylabel('w*100')
title('Phase-Plane Plot (MLE)')

[V,w]=meshgrid(-84:10:120,0:0.05:1);     %For quiver plot
dV=(1/C)*(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Eca)-Gk*w.*(V-Ek)-Gl*(V-El));
dw=phi*(0.5*(1+tanh((V-V3)/V4))-w).*cosh((V-V3)/(2*V4));
quiver(V,w*100,dV,dw*100);
%}
%
%% 3. Jacobian and Eigenvalues
syms V w;
jac=jacobian([(1/C)*(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Eca)-Gk*w.*(V-Ek)-Gl*(V-El)), phi*(0.5*(1+tanh((V-V3)/V4))-w).*cosh((V-V3)/(2*V4))], [V,w]);
disp(['Iext = ' num2str(Iext) '     Veq = ' num2str(double(Veq)) '     weq = ' num2str(double(weq))]);
jac=double(subs(jac,[V,w],[Veq1,weq1]));
eigenvalues=eig(jac)
%}
%
%% 5. Action Potential versus phi
figure(2)
hold on
fplot(@(V) V_nc(V),[-84 120 0 1],'k');
fplot(@(V) w_nc(V),[-84 120 0 1],'k');
xlabel('V(mV)')
ylabel('w')
title('Phase-Plane Plot (MLE)')

phi_arr=[0.01 0.02 0.04];
for i=1:3
    tspan=[0 300];
    x0=[-13.9;weq];
    phi=phi_arr(i);
    [t,x]=ode15s(@MLE,tspan,x0);
    plot(x(:,1),x(:,2));
    figure(3)
    hold on
    plot(t,x(:,1));                             % 2nd method to find eqm pt - see trajectory
    figure(2)
end

h=get(gca,'Children');
legend([h(3) h(2) h(1)],['\phi=' num2str(phi_arr(1))],['\phi=' num2str(phi_arr(2))],['\phi=' num2str(phi_arr(3))])
figure(3)
xlabel('t (ms)')
ylabel('V(t)')
title('Action Potential for MLE')
legend('\phi=0.01','\phi=0.02','\phi=0.04')
phi=0.02;
%}
%
%% 6. Action potential of depolarizing current pulses
figure(4)
hold on
fplot(@(V) V_nc(V),[-84 120 0 1],'k');
fplot(@(V) w_nc(V),[-84 120 0 1],'k');
xlabel('V(mV)')
ylabel('w')
title('Phase-Plane Plot (MLE)')
V0=[-15 -14.8];
tspan=[0 300];
for i=1:2
    x0=[V0(i);weq];
    [~,x]=ode15s(@MLE,tspan,x0);
    plot(x(:,1),x(:,2));
end
h=get(gca,'Children');
legend([h(2) h(1)],'V_0=-15mV','V_0=-14.8mV')

V0=[-25:5:-15 -14.99:0.01:-13 -12.9:0.1:-5];
Vpp=zeros(1,length(V0));

for i=1:length(V0)
    [~,x]=ode15s(@MLE,tspan,[V0(i);weq]);
    Vpp(i)=max(x(:,1));
end
figure(5)
plot(V0,Vpp);
ylabel('Amplitude')
xlabel('Initial Voltage (V_0)')
title('Action potential amplitude (MLE)')
%}
%
%% 7. Effect of External current
Iext=86;
syms V w;
V_nc=@(V) (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gl*(V-El))/(Gk*(V-Ek));
w_nc=@(V) (0.5*(1+tanh((V-V3)/V4)));
figure(6)
hold on
fplot(@(V) V_nc(V),[-84 120 0 1],'k');
fplot(@(V) w_nc(V),[-84 120 0 1],'k');
xlabel('V(mV)')
ylabel('w')
title('Phase-Plane Plot (MLE)')

%Three initial conditions
V0=zeros(1,3);
w0=zeros(1,3);
V0(1)=Veq;w0(1)=weq;
V_nc=(Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gk*w*(V-Ek)-Gl*(V-El))==0;
w_nc=(0.5*(1+tanh((V-V3)/V4))-w)==0;
[Veq,weq]=solve([V_nc,w_nc],[V,w]);
V0(2)=double(Veq);w0(2)=double(weq);
V0(3)=-27.9;w0(3)=0.17;

tspan=[0 300];
for i=1:3
    x0=[V0(i);w0(i)];
    [t,x]=ode15s(@MLE,tspan,x0);
    plot(x(:,1),x(:,2));
    figure(7)
    hold on
    plot(t,x(:,1));                             % 2nd method to find eqm pt - see trajectory
    figure(6)
end
h=get(gca,'Children');
legend([h(3) h(2) h(1)],'X_0=Old equilibrium','X_0=New equilibrium','X_0=(-27.9,0.17)')
figure(7)
h=get(gca,'Children');
legend([h(3) h(2) h(1)],'X_0=Old equilibrium','X_0=New equilibrium','X_0=(-27.9,0.17)')
title('Action Potentials (I_e_x_t=86\muA/cm^2)')

syms V w;
jac=jacobian([(1/C)*(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Eca)-Gk*w.*(V-Ek)-Gl*(V-El)), phi*(0.5*(1+tanh((V-V3)/V4))-w).*cosh((V-V3)/(2*V4))], [V,w]);
jac=double(subs(jac,[V,w],[Veq,weq]));
disp(['Iext = ' num2str(Iext) '     Veq = ' num2str(double(Veq)) '     weq = ' num2str(double(weq))]);
eigenvalues=eig(jac)
%}
%
%% 8. Finding the UPO to characterize bistability
V_nc=@(V) (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gl*(V-El))/(Gk*(V-Ek));
w_nc=@(V) (0.5*(1+tanh((V-V3)/V4)));
for i=1:2
    figure(10-i);
    hold on
    fplot(@(V) V_nc(V),[-84 120 0 1],'k');
    fplot(@(V) w_nc(V),[-84 120 0 1],'k');
    xlabel('V(mV)')
    ylabel('w')
    title('Phase-Plane Plot (MLE)')
end

tspan1=[3000 0];                        % Run backwards in time from inf to find limiting UPO
x0=[double(Veq);double(weq)+0.01];                 % Small displacement from equilibrium point
[~,x]=ode15s(@MLE,tspan1,x0);
plot(x(:,1),x(:,2));
h=get(gca,'Children');
legend(h(1),'Reverse-time trajectory for UPO')

figure(9)
K_ms = convhull(x(:,1),x(:,2));        % Take convex hull to get UPO
plot(x(K_ms,1),x(K_ms,2));
x0=[-17.2;0.16];
[~,x]=ode15s(@MLE,tspan,x0);
plot(x(:,1),x(:,2));
x0=[-16.9;0.16];
[~,x]=ode15s(@MLE,tspan,x0);
plot(x(:,1),x(:,2));
title('Bistability in MLE')
h=get(gca,'Children');
legend([h(3) h(2) h(1)],'UPO','Stable Spiral','Limit Cycle (SPO)');

figure(10)
V0=-22.5:0.1:-17.5;
Vpp=zeros(1,length(V0));
for i=1:length(V0)
    [~,x]=ode15s(@MLE,[0 1000],[V0(i);0.1376]);
    Vpp(i)=max(x(:,1));
end
plot(V0,Vpp);
ylabel('Amplitude')
xlabel('Initial Voltage (V_0)')
title('Action potential amplitude - Bistability(MLE)')
%}
%
%% 9. Effect of external current on rate of firing of action potentials
I_arr=[80 86 90];
figure(11)
tspan=[0 1000];
for i=1:3
    Iext=I_arr(i);
    syms V w;
    V_nc=@(V) (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gl*(V-El))/(Gk*(V-Ek));
    w_nc=@(V) (0.5*(1+tanh((V-V3)/V4)));
    subplot(2,2,i);
    hold on
    fplot(@(V) V_nc(V),[-84 120 0 1],'k');
    fplot(@(V) w_nc(V),[-84 120 0 1],'k');
    xlabel('V(mV)')
    ylabel('w')
    title('Phase-Plane Plot (MLE)')

    V_nc=(Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gk*w*(V-Ek)-Gl*(V-El))==0;
    w_nc=(0.5*(1+tanh((V-V3)/V4))-w)==0;
    [Veq,weq]=solve([V_nc,w_nc],[V,w]);
    
    jac=jacobian([(1/C)*(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Eca)-Gk*w.*(V-Ek)-Gl*(V-El)), phi*(0.5*(1+tanh((V-V3)/V4))-w).*cosh((V-V3)/(2*V4))], [V,w]);
    jac=double(subs(jac,[V,w],[Veq,weq]));
    disp(['Iext = ' num2str(Iext) '     Veq = ' num2str(double(Veq)) '     weq = ' num2str(double(weq))]);
    eigenvalues=eig(jac)
    
    x0=[double(Veq);double(weq)+0.01];
    [~,x]=ode15s(@MLE,tspan,x0);
    plot(x(:,1),x(:,2),'r');
    h=get(gca,'Children');
    legend(h(1),['Trajectory for I_e_x_t=' num2str(I_arr(i)) '\muA/cm^2'])
end

I_arr=[80:4:88 88.1:0.1:91.9 92:4:100];
rate_firing=zeros(1,length(I_arr));

tspan=[1000 2000];
for i=1:length(I_arr)
    Iext=I_arr(i);
    V_nc=(Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gk*w*(V-Ek)-Gl*(V-El))==0;
    w_nc=(0.5*(1+tanh((V-V3)/V4))-w)==0;
    [Veq,weq]=solve([V_nc,w_nc],[V,w]); x0=[double(Veq);double(weq)+0.01];
    [~,x]=ode15s(@MLE,tspan,x0);
    rate_firing(i)=sum(x(1:length(x(:,1))-1) < 0 & x(2:length(x(:,1))) > 0 );
end
figure(12)
plot(I_arr,rate_firing)
xlabel('Applied current (\muA/cm^2)')
ylabel('Rate of firing (Hz)')
title('Rate of firing of action potentials (MLE 1^s^t set of parameters)')
%}
%
%% 10. Second set of parameters
Gca=4; Gk=8.0; Gl=2; Eca=120; Ek=-84; El=-60;             %values of part 10
phi=0.0667; V1=-1.2; V2=18; V3=12; V4=17.4; C=20; Iext=30;
V_nc=@(V) (Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gl*(V-El))/(Gk*(V-Ek));  % for nullclines
w_nc=@(V) (0.5*(1+tanh((V-V3)/V4)));
figure(13)
hold on
fplot(@(V) V_nc(V),[-84 120 0 1],'k');
fplot(@(V) w_nc(V),[-84 120 0 1],'k');
xlabel('V(mV)')
ylabel('w')
title('Phase-Plane Plot (MLE 2^n^d set of parameters)')
%
syms V w;
V_nc=(Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gk*w*(V-Ek)-Gl*(V-El))==0;
w_nc=(0.5*(1+tanh((V-V3)/V4))-w)==0;

Veq=sym('Veq',[1 20]);
weq=sym('weq',[1 20]);
for i=1:length(Veq)
    [Veq(i),weq(i)]=vpasolve([V_nc,w_nc],[V,w],'random',true);
end
eq=unique([double(Veq'),double(weq')],'rows');

disp('For the 2nd set of MLE parameters:');
for i=1:size(eq,1)
    jac=jacobian([(1/C)*(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Eca)-Gk*w.*(V-Ek)-Gl*(V-El)), phi*(0.5*(1+tanh((V-V3)/V4))-w).*cosh((V-V3)/(2*V4))], [V,w]);
    jac=double(subs(jac,[V,w],[eq(i,1),eq(i,2)]));
    disp(['Veq = ' num2str(eq(i,1)) '     weq = ' num2str(eq(i,2))]);
    eigenvalues=eig(jac)
    h=plot(Veq,weq,'y+');
    set(h,'linewidth',2)
    text(eq(i,1),eq(i,2),['  \leftarrow(' num2str(round(eq(i,1),3)) ',' num2str(round(eq(i,2),3)) ')'])
end

% 3eqm points -> Stable, Saddle node, Unstable Spiral
%Unstable manifolds
tspan=[0 2000];
x0=[eq(2,1)-0.01;eq(2,2)+0.001];         %Saddle node
[~,x]=ode15s(@MLE,tspan,x0);
plot(x(:,1),x(:,2),'r');

x0=[eq(2,1)+0.01;eq(2,2)-0.001];         %Saddle node
[~,x]=ode15s(@MLE,tspan,x0);
plot(x(:,1),x(:,2),'r');

x0=[eq(3,1)+0.01;eq(3,2)];         %Unstable spiral
[~,x]=ode15s(@MLE,tspan,x0);
plot(x(:,1),x(:,2),'r');

% Remaining stable manifolds
tspan=[2000 0];
x0=[eq(2,1)+0.1;eq(2,2)+0.001];         %Saddle node
[~,x]=ode15s(@MLE,tspan,x0);
plot(x(:,1),x(:,2),'b');

tspan=[100 0];
x0=[eq(2,1);eq(2,2)-0.001];         %Saddle node
[~,x]=ode15s(@MLE,tspan,x0);
plot(x(:,1),x(:,2),'b');

h=get(gca,'Children');
legend([h(1) h(3)],'Stable manifold', 'Unstable manifold')
annotation('textbox',[0.75 0.7 0.11 0.1],'String','Eqm. points (left to right): Stable, Saddle Node, Unstable Spiral','FitBoxToText','off')
%}
%
%% 11. Change of equilibrium points and rate of firing with current
I_arr=[30 35:2:39 39.1:.1:41 42:4:50];
eqm=cell(length(I_arr),3);
syms V w;
n=20;
rate_firing=zeros(1,length(I_arr));
tspan=[1000 2200];

disp('Terminology: Stable->S  Unstable->US  Saddle Node->SN  Stable Spiral->SS  Unstable Spiral->US');
for i=1:length(I_arr)
    Iext=I_arr(i);
    V_nc=(Iext-Gca*0.5*(1+tanh((V-V1)/V2))*(V-Eca)-Gk*w*(V-Ek)-Gl*(V-El))==0;
    w_nc=(0.5*(1+tanh((V-V3)/V4))-w)==0;
    Veq=sym('Veq',[1 n]);
    weq=sym('weq',[1 n]);
    for j=1:length(Veq)
        [Veq(j),weq(j)]=vpasolve([V_nc,w_nc],[V,w],'random',true);
    end
    eq=unique([double(Veq'),double(weq')],'rows');
    if(size(eq,1)==1)
        n=1;
    end
    
    for j=1:size(eq,1)
        jac=jacobian([(1/C)*(Iext-Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Eca)-Gk*w.*(V-Ek)-Gl*(V-El)), phi*(0.5*(1+tanh((V-V3)/V4))-w).*cosh((V-V3)/(2*V4))], [V,w]);
        jac=double(subs(jac,[V,w],[eq(j,1),eq(j,2)]));
        t=eig(jac);
        if(all(imag(t)==0)&&all(real(t)>0))
            eqm(i,j+(3-size(eq,1)))={'U'};
        elseif(all(imag(t)==0)&&all(real(t)<0))
            eqm(i,j+(3-size(eq,1)))={'S'};
        elseif(all(imag(t)==0)&&(real(t(1))*real(t(2))<0))
            eqm(i,j+(3-size(eq,1)))={'SN'};
        elseif(all(imag(t)~=0)&&all(real(t)>0))
            eqm(i,j+(3-size(eq,1)))={'US'};
        elseif(all(imag(t)~=0)&&all(real(t)<0))
            eqm(i,j+(3-size(eq,1)))={'SS'};
        end
    end
    
    x0=[-19.5;0.026];
    [t,x]=ode15s(@MLE,tspan,x0);
    rate_firing(i)=sum(x(1:length(x(:,1))-1) < 0 & x(2:length(x(:,1))) > 0 );
end

eqm1=eqm(:,1);eqm2=eqm(:,2);eqm3=eqm(:,3);
EqmVsIext=table(eqm1,eqm2,eqm3,'RowNames',arrayfun(@num2str, I_arr, 'Uniform', false)')

figure(14)
plot(I_arr,rate_firing);
xlabel('Applied current (\muA/cm^2)')
ylabel('Rate of firing (Hz)')
title('Rate of firing of action potentials (MLE 2^n^d set of parameters)')
%}

end


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