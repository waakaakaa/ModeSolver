function [eeps,nx,ny,mStEnd] = GenerateEps(vT,vN,vW,ncld,wid,dx,dy)%ȡ�����ҶԳƵĶ�㲨���ṹ�������ʷֲ�����
%vT,vN,vW�ֱ���ÿһ���ȣ������ʣ��Ϳ��
%wid���ܿ�ȣ�vW������wid,�ᱻ�Զ��ضϣ���<0����Ϊwid
%mStEnd 4Ԫ���У���ʾeeps��ÿһ�����ʼ����ֹ���
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
    mStEnd(1,i) = round((nx-nxi)/2)+1;%x������ʼ���
    mStEnd(2,i) = round((nx+nxi)/2);%x������ֹ���
    mStEnd(3,i) = n1;%y������ʼ���
    mStEnd(4,i) = n1+nyi-1;%y������ֹ���
    n1 = n1+nyi;
end