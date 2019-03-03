function dz = double_gyre(t,z)

n = numel(z)/2;
x = z(1:n,1);
y = z(n+1:2*n,1);

st = ((t>0)&(t<1)).*t.^2.*(3-2*t) + (t>1)*1;

dxPsi_P = 2*pi*cos(2*pi*x).*sin(pi*y);
dyPsi_P = pi*sin(2*pi*x).*cos(pi*y);
dxPsi_F = pi*cos(pi*x).*sin(2*pi*y);
dyPsi_F = 2*pi*sin(pi*x).*cos(2*pi*y);
dxPsi = (1-st).*dxPsi_P + st.*dxPsi_F;
dyPsi = (1-st).*dyPsi_P + st.*dyPsi_F;

dz(1:n,1) = -dyPsi;
dz(n+1:2*n,1) = dxPsi;



