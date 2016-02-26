function [vEpsl] = GetEpslProf1D(vThn,vN,dSlot,dx,varargin)%epsl profile in waveguide section with slot 
%dslot is the depth of the slot
%
%                  vN(1) vThn(1)  
%     ------------    -----------------------------
%                 |  |         
%                 |  |   vN(2)   vThn(2)     
%     ____________|  |_____________________________ 
%                 |__|   vN(3)   vThn(3)  
%     _____________________________________________ 
%                        vN(4)   vThn(4)  
%     _____________________________________________ 
% 
%                        .
%                        .
%                        .
%     _____________________________________________  
%                   vN(n)       vThn(n)  
%     _____________________________________________  

if(nargin>4)
    iStart = varargin{1};%�̲۵���ʼ��
else
    iStart = 2;%����Ĭ�ϴӵڶ��㿪ʼ��(��һ��Ϊ���ʲ�)
end
n = length(vThn);
vn = round(vThn/dx);
vEpsl = [];
if(dSlot ~= Inf)   % slot�ĺ�ȣ������ǿ��
    k = 0;
    for i = 1:n;
        vEpsl = [vEpsl,vN(i)^2*ones(1,vn(i))];
        if(i<iStart)
            k = k+vn(i);
        end
    end;
    nslot = round(dSlot/dx);
    vEpsl(k+1:k+nslot) = vN(1)^2*ones(1,nslot);
else
    vEpsl = vN(1)^2*ones(1,sum(vn));
end;
