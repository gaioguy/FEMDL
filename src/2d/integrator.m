function [xt,yt] = integrator(v,x0,y0,tspan)

options = odeset('RelTol',1e-3,'AbsTol',1e-3);

[~,F] = ode45(v,tspan,[x0(:);y0(:)],options);

if numel(tspan)==2,
    xt = F([1,end],1:end/2);
    yt = F([1,end],end/2+1:end);
else
    xt = F(:,1:end/2);
    yt = F(:,end/2+1:end);
end

