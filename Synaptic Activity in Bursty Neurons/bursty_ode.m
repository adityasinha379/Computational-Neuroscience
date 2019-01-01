function dx = bursty_ode(x,Iext,N)
V=x(1);

xinf=0.5*(1+tanh(N.ko.*(V-N.Vth)))';      % Ca(n) K(m) KS(c)
tau=sech(N.ko(2)/2.*(V-N.Vth(2:3)))';         % K(m) KS(c)
% Update outputs
dx(1)=1/N.C*(Iext-sum(N.g.*[xinf(1) x(2:3)' 1].*(V-N.E)));
dx(2:3)=[N.eps; N.del]./tau.*(xinf(2:3)-x(2:3));
end

