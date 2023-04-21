t1=0; 
t2=1;
h=0.01;
n = ceil((t2-t1)/h);
tt = zeros(n);
t = zeros(n);
u0=0; 
z0 = zeros(n);
u = zeros(1,n);
z0(1)=0.3;
z0(2)=0.31;

thresh = 1e-5;
  
     h = (t2-t1)/n;     
     x = t1:h:t2;


[t,tt] = RK4(t1, u0, z0(1), h, t2);
u(1) = tt(end);
for j=1:5
    h=h/2;
 i = 1;
while true
i = i+ 1;
[t,tt] = RK4(t1, u0, z0(i), h, t2);
u(i) = tt(end);
 

if abs(1-u(i))< thresh     
break
end

z0(i+1) = z0(i) - ((u(i)-1)*(z0(i-1) - z0(i)))/((u(i-1)-1)-(u(i)-1));
hold on
plot(t,tt)
hold off
end
end
