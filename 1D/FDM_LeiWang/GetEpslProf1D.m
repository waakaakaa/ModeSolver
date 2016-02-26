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
    iStart = varargin{1};%刻槽的起始层
else
    iStart = 2;%否则默认从第二层开始刻(第一层为介质层)
end
n = length(vThn);
vn = round(vThn/dx);
vEpsl = [];
if(dSlot ~= Inf)   % slot的厚度，而不是宽度
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
