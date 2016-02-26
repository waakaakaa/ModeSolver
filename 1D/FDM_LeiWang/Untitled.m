n=1.4591;
dn=1.4582-n;
na = 1.4565;
dna=1.4554-na;
lamda=1.49;
dlam=0.06;
m = 10;
A=[n,na;n+dn,na+dna];
B = [m*lamda;m*(lamda+dlam)];
X = inv(A)*B


n=2.978936;
dn=2.954395-n;
na = 2.554248;
dna=2.494389-na;
lamda=1.285;
dlam=0.05;
m = 10;
A=[n,na;n+dn,na+dna];
B = [m*lamda;m*(lamda+dlam)];
X = inv(A)*B
