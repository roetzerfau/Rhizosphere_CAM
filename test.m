tspan = [0 5];
y0 = [1 1];
[t,y] = ode45(@(t,y) odefcn(t,y), tspan, y0);
plot(t,y,'-o')
legend()

function dydt = odefcn(t,y)
if(y(1)<0)
    y(1)
end
  dydt = zeros(2,1);
  dydt(1) = - y(1);
  dydt(2) = - 3*y(2);
end