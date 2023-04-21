function [t,u] = RK4(t1,u0,z0,h,t2)
n = ceil((t2-t1)/h);  % number of elements
t = zeros(1,n);   % defining the t array
u = zeros(1,n);   % defining the y array
z = zeros(n,1);   % defining the z array
% defining all the coefficients for RK3 method
a1=1/6;
a2=2/6;
a3=2/6;
a4=1/6;
p1=1/2;
p2=1/2;
p3=1;
q11=1/2;
q21=0;
q22=1/2;
q31=0;
q32=0;
q33=1;

% defining the initial values of t,y,z and final value of t
t(1) = t1; 
u(1) = u0;
z(1) = z0;
t(n) = t2;
h = (t2-t1)/n;  % step size such that the final value coincides with t2

for i = 1:n  % this loop calculates t,y,z arrays using rk3 method

    k1 = dudt(t(i),u(i),z(i));
    l1 = dzdt(t(i),u(i),z(i));
    k2 = dudt(t(i)+p1*h,u(i)+k1*q11*h,z(i)+l1*q11*h);
    l2 = dzdt(t(i)+p1*h,u(i)+k1*q11*h,z(i)+l1*q11*h);
    k3 = dudt(t(i)+p2*h,u(i)+k1*q21*h+k2*q22*h,z(i)+l1*q21*h+l2*q22*h);
    l3 = dzdt(t(i)+p2*h,u(i)+k1*q21*h+k2*q22*h,z(i)+l1*q21*h+l2*q22*h);
    k4 = dudt(t(i)+p3*h,u(i)+k1*q31*h+k2*q32*h+k3*q33*h,z(i)+l1*q31*h+l2*q32*h+l3*q33*h);
    l4 = dudt(t(i)+p3*h,u(i)+k1*q31*h+k2*q32*h+k3*q33*h,z(i)+l1*q31*h+l2*q32*h+l3*q33*h);
    phi1 = a1*k1+a2*k2+a3*k3+a4*k4;
    phi2 = a1*l1+a2*l2+a3*l3+a4*l4;
    

    t(i+1) = t(i) + h;
    u(i+1) = u(i) + h*phi1;
    z(i+1) = z(i) + h*phi2;
    
end
return
end