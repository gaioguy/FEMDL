function [xt,yt,zt] = integrator3(v,x0,y0,z0,tspan)

options = odeset('RelTol',1e-3,'AbsTol',1e-3);

[T,F] = ode45(v,tspan,[x0(:);y0(:);z0(:)],options);

F = mod(F,1);
if numel(tspan)==2
    xt = F([1,end],1:end/3);
    yt = F([1,end],end/3+1:end*2/3);
    zt = F([1,end],end*2/3+1:end);
else
    xt = F(:,1:end/3);
    yt = F(:,end/3+1:end*2/3);
    zt = F(:,end*2/3+1:end);
end

