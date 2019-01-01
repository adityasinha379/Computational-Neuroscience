function tk = PN(t,v,b,kappa,delta)

N=size(v,1);
tk = cell(N,1);
dt = diff(t(1:2));

for i=1:N
    temp = 0;
    for j = 1:length(t)
        temp = temp + dt/kappa(i)*(v(i,j)+b(i));
        if temp >= delta(i)
            tk{i} = [tk{i} t(j)];
            temp = temp - delta(i);
        end
    end
end

end

