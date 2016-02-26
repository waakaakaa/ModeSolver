function Gama = GetConfineBylayerNumbers(mE,nx,ny,iStart,iEnd,mStEnd)%获得iStart layer 到iEnd layer之间光能量占总能量的比例
%mStEnd 是每一层起始和终止格点的位置
i1 = mStEnd(3,iStart);
i2 = mStEnd(4,iEnd);
Gama = sum(sum(abs(mE(:,i1:i2).^2)))/sum(sum(abs(mE(:,:).^2)));
