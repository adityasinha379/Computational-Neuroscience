function tk = OSN(t,u,h,b,kappa,delta)

dt = diff(t(1:2));
M=size(h,1);
tk = cell(M,1);


for i=1:M
    v = conv(u,h(i,:))*dt;
    v = v(1:length(t));
    temp = 0;
    for j = 1:length(t)
        temp = temp + dt/kappa(i)*(v(j)+b(i));
        if temp >= delta(i)
            tk{i} = [tk{i} t(j)];
            temp = temp - delta(i);
        end
    end
end

end

