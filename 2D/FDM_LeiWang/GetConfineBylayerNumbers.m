function Gama = GetConfineBylayerNumbers(mE,nx,ny,iStart,iEnd,mStEnd)%���iStart layer ��iEnd layer֮�������ռ�������ı���
%mStEnd ��ÿһ����ʼ����ֹ����λ��
i1 = mStEnd(3,iStart);
i2 = mStEnd(4,iEnd);
Gama = sum(sum(abs(mE(:,i1:i2).^2)))/sum(sum(abs(mE(:,:).^2)));
