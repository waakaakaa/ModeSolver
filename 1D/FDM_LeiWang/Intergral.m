function y = Intergral(x1,y1,x2,y2,dx)%abs^2(Quad(y1*conj*y2))/(Quad(abs^2(y2))*Quad(abs^2(y1))
min1 = min(x1);
min2 = min(x2);
max1 = max(x1);
max2 = max(x2);

n0 = round(abs((min2-min1)/dx));
if(min1 < min2)
    y2 = [zeros(1,n0),y2];
elseif(min1 > min2)
    y1 = [zeros(1,n0),y1];
end;

n0 = round(abs((max2-max1)/dx));
if(max1 > max2)
    y2 = [y2,zeros(1,n0)];
elseif(max1 < max2)
    y1 = [y1,zeros(1,n0)];
end;
n1 = length(y1);
n2 = length(y2);
if(n1>n2)
    y1 = y1(1:n2);
elseif(n2>n1)
     y2 = y2(1:n2);
end;
     y = abs(sum(y1.*conj(y2)))^2/(abs(sum(y2.*conj(y2)))*abs(sum(y1.*conj(y1))));