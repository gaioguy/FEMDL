function dy = bickleyjet_vf(t,y)

U=62.66e-6;
L=1770e-3;
c2=.205*U;
c3=.461*U;
eps1=0.0075;
eps2=0.15;
eps3=0.3;
r0=6371e-3;
k1=2/r0;
k2=4/r0;
k3=6/r0;
c1=c3+((sqrt(5)-1)/2)*(k2/k1)*(c2-c3);

N = round(length(y)/2);
dy = zeros(2*N,1);

dy(1:N,1) = U*sech(y(N+1:2*N,1)/L).^2+(2*eps1*U*cos(k1*(y(1:N,1)-c1*t))+...
    2*eps2*U*cos(k2*(y(1:N,1)-c2*t))+...
    2*eps3*U*cos(k3*(y(1:N,1)-c3*t))).*tanh(y(N+1:2*N,1)/L).*sech(y(N+1:2*N,1)/L).^2;

dy(N+1:2*N,1) = -(eps1*k1*U*L*sin(k1*(y(1:N,1)-c1*t))+...
    eps2*k2*U*L*sin(k2*(y(1:N,1)-c2*t))+...
    eps3*k3*U*L*sin(k3*(y(1:N,1)-c3*t))).*sech(y(N+1:2*N,1)/L).^2;

