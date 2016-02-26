function [eeps,nx,ny] = getEpsl(vThickness,vNeffact,vWidth,totalwid,dx,dy)%取得左右对称的多层波导结构的折射率分布矩阵

n = length(vThickness);
vnxi =  round(vWidth/dx);
nx = round(totalwid/dx);
vnyi = round(vThickness/dy);
ny = sum(vnyi);
eeps = ones(nx,ny);
n1 = 1;
for i=1:n
    nxi = vnxi(i);
    nyi = vnyi(i);
    eeps(round((nx-nxi)/2)+1:round((nx+nxi)/2),n1:n1+nyi-1) = vNeffact(i).^2*ones(nxi,nyi);
    n1 = n1 + nyi;
end