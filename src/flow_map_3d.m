function y = flow_map3(v,x,tspan)

[Fx,Fy,Fz] = integrator_3d(v,x(:,1),x(:,2),x(:,3),tspan);
y = [Fx' Fy' Fz']; 

end