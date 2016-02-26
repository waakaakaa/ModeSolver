function [vEpsl] = getEpslFromThickAndN(vThickness,vNeffAct,dx)
n = length(vThickness);
vn = round(vThickness/dx);
vEpsl = [];
for i = 1:n
    vEpsl = [vEpsl, vNeffAct(i)^2 * ones(1,vn(i))];
end
