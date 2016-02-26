function [eeps,nx,ny,mStEnd] = GenerateEps(vT,vN,vW,ncld,wid,dx,dy)%取得左右对称的多层波导结构的折射率分布矩阵
%vT,vN,vW分别是每一层厚度，折射率，和宽度
%wid是总宽度，vW若大于wid,会被自动截断，若<0，则为wid
%mStEnd 4元序列，表示eeps中每一层的起始和终止序号
n = length(vT);
vnxi =  round(vW/dx);
nx = round(wid/dx);
vnyi = round(vT/dy);
ny = sum(vnyi);
eeps = ncld.^2*ones(nx,ny);
n1 = 1;
mStEnd = zeros(4,n);
for i=1:n
    if(vnxi(i)<0 || vnxi(i)>nx)
        nxi = nx;
    else
        nxi = vnxi(i);
    end
    nyi = vnyi(i);
    eeps(round((nx-nxi)/2)+1:round((nx+nxi)/2),n1:n1+nyi-1) = vN(i).^2*ones(nxi,nyi);
    mStEnd(1,i) = round((nx-nxi)/2)+1;%x方向起始序号
    mStEnd(2,i) = round((nx+nxi)/2);%x方向终止序号
    mStEnd(3,i) = n1;%y方向起始序号
    mStEnd(4,i) = n1+nyi-1;%y方向终止序号
    n1 = n1+nyi;
end