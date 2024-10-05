clear all; close all; clc; 
%calculate the free energy of the new solid solution for binary system of
%FCC structure

%binary system
%tetrahedron

%% 0. constants
global kB N T x 
kB=1.380649*10^-23;  %boltzman constant, unit:J/K
N=6.02*10^23;        %number of particles, the size of the system, here we use 1 mole

T=1:0.001:2.5;
T=T';
%% 1. model setting


%pair energy eV
E11=0;
E00=0;
E10=0;
E01=0;

%F(T,x)

%x=0.49;
parameter1=zeros(1501,1); %xxx
parameter2=zeros(1501,1); %X

minF=zeros(1501,1);
minF=minF+100;
%xxx=-0.13;
i=0;

x=0.5;
tic
i=i+1;

for xxx=-7:0.0001:0
%for yyy=-5:0.1:0


L1X=xxx;
L2X=xxx;
L3X=-xxx;
L4X=-xxx;


lambda=0;


%constants/parameters
%cluster energy, unit?

E0000=0;
E1111=0;
E1110=-3*kB;
E1101=E1110;
E1100=-4*kB;
E1011=E1110;
E1010=E1100;
E1001=E1100;
E1000=-3*kB;
E0111=E1110;
E0110=E1100;
E0101=E1100;
E0100=E1000;
E0011=E1100;
E0010=E1000;
E0001=E1000;

%must check the order-disorder transition's property, continuous, but derivative has the jump.
%use Landau's chapter 14 of volume 5 to check it.

%here Eorder must be a ordered phase energy. However, if the composition is
%not exactly at some ratio, it must 



vv(1,1)=L1X;
vv(1,2)=L2X;
vv(1,3)=L3X;
vv(1,4)=L4X;
vv(1,5)=lambda;



part4=exp(-E0000/kB./T)+exp(L1X-E1000/kB./T)+exp(L2X-E0100/kB./T)+exp(L3X-E0010/kB./T)+exp(L4X-E0001/kB./T)+ ...,
    exp(L1X+L2X-E1100/kB./T)+exp(L2X+L3X-E0110/kB./T)+exp(L3X+L4X-E0011/kB./T)+exp(L1X+L4X-E1001/kB./T)+exp(L1X+L3X-E1010/kB./T)+exp(L2X+L4X-E0101/kB./T)+ ...,
    exp(L1X+L2X+L3X-E1110/kB./T)+exp(L2X+L3X+L4X-E0111/kB./T)+exp(L1X+L3X+L4X-E1011/kB./T)+exp(L1X+L2X+L4X-E1101/kB./T)+ ...,
    exp(L1X+L2X+L3X+L4X-E1111/kB./T);

X1=(exp(L1X-E1000/kB./T)+exp(L1X+L2X-E1100/kB./T)+exp(L1X+L4X-E1001/kB./T)+exp(L1X+L3X-E1010/kB./T)+ ...,
    exp(L1X+L2X+L3X-E1110/kB./T)+exp(L1X+L3X+L4X-E1011/kB./T)+exp(L1X+L2X+L4X-E1101/kB./T)+exp(L1X+L2X+L3X+L4X-E1111/kB./T))./part4;

X2=(exp(L2X-E0100/kB./T)+exp(L1X+L2X-E1100/kB./T)+exp(L2X+L3X-E0110/kB./T)+exp(L2X+L4X-E0101/kB./T)+ ...,
    exp(L1X+L2X+L3X-E1110/kB./T)+exp(L2X+L3X+L4X-E0111/kB./T)+exp(L1X+L2X+L4X-E1101/kB./T)+exp(L1X+L2X+L3X+L4X-E1111/kB./T))./part4;

X3=(exp(L3X-E0010/kB./T)+exp(L2X+L3X-E0110/kB./T)+exp(L3X+L4X-E0011/kB./T)+exp(L1X+L3X-E1010/kB./T)+ ...,
    exp(L1X+L2X+L3X-E1110/kB./T)+exp(L2X+L3X+L4X-E0111/kB./T)+exp(L1X+L3X+L4X-E1011/kB./T)+exp(L1X+L2X+L3X+L4X-E1111/kB./T))./part4;

X4=(exp(L4X-E0001/kB./T)+exp(L3X+L4X-E0011/kB./T)+exp(L1X+L4X-E1001/kB./T)+exp(L2X+L4X-E0101/kB./T)+ ...,
    exp(L2X+L3X+L4X-E0111/kB./T)+exp(L1X+L3X+L4X-E1011/kB./T)+exp(L1X+L2X+L4X-E1101/kB./T)+exp(L1X+L2X+L3X+L4X-E1111/kB./T))./part4;
	
	
	
	
	
	
X1111=exp(L1X+L2X+L3X+L4X-E1111/kB./T)./part4;
X1110=exp(L1X+L2X+L3X-E1110/kB./T)./part4;
X1101=exp(L1X+L2X+L4X-E1101/kB./T)./part4;
X1011=exp(L1X+L3X+L4X-E1011/kB./T)./part4;
X0111=exp(L2X+L3X+L4X-E0111/kB./T)./part4;
X1100=exp(L1X+L2X-E1100/kB./T)./part4;
X0110=exp(L2X+L3X-E0110/kB./T)./part4;
X0101=exp(L2X+L4X-E0101/kB./T)./part4;
X1010=exp(L1X+L3X-E1010/kB./T)./part4;
X1001=exp(L1X+L4X-E1001/kB./T)./part4;
X0011=exp(L3X+L4X-E0011/kB./T)./part4;
X1000=exp(L1X-E1000/kB./T)./part4;
X0100=exp(L2X-E0100/kB./T)./part4;
X0010=exp(L3X-E0010/kB./T)./part4;
X0001=exp(L4X-E0001/kB./T)./part4;
X0000=exp(-E0000/kB./T)./part4;








%A=12,B=13,C=14,D=23,E=24,F=34

X00A=X0000+X0001+X0010+X0011;
X00B=X0000+X0001+X0100+X0101;
X00C=X0000+X0010+X0100+X0110;
X00D=X0000+X0001+X1000+X1001;
X00E=X0000+X0010+X1000+X1010;
X00F=X0000+X1000+X0100+X1100;


X11A=X1100+X1101+X1110+X1111;
X11B=X1010+X1011+X1110+X1111;
X11C=X1001+X1011+X1101+X1111;
X11D=X0110+X0111+X1110+X1111;
X11E=X0101+X0111+X1101+X1111;
X11F=X0011+X1011+X0111+X1111;


X10A=X1000+X1001+X1010+X1011;
X10B=X1000+X1001+X1100+X1101;
X10C=X1000+X1010+X1100+X1110;
X10D=X0100+X0101+X1100+X1101;
X10E=X0100+X0110+X1100+X1110;
X10F=X0010+X1010+X0110+X1110;


X01A=X0100+X0101+X0110+X0111;
X01B=X0010+X0011+X0110+X0111;
X01C=X0001+X0011+X0101+X0111;
X01D=X0010+X0011+X1010+X1011;
X01E=X0001+X0011+X1001+X1011;
X01F=X0001+X1001+X0101+X1101;



FF2=(X00A.*(E00+kB.*T.*log(X00A))+X11A.*(E11+kB.*T.*log(X11A))+X10A.*(E10+kB.*T.*log(X10A))+X01A.*(E01+kB.*T.*log(X01A)))+ ...,
    (X00B.*(E00+kB.*T.*log(X00B))+X11B.*(E11+kB.*T.*log(X11B))+X10B.*(E10+kB.*T.*log(X10B))+X01B.*(E01+kB.*T.*log(X01B)))+ ...,
    (X00C.*(E00+kB.*T.*log(X00C))+X11C.*(E11+kB.*T.*log(X11C))+X10C.*(E10+kB.*T.*log(X10C))+X01C.*(E01+kB.*T.*log(X01C)))+ ...,
    (X00D.*(E00+kB.*T.*log(X00D))+X11D.*(E11+kB.*T.*log(X11D))+X10D.*(E10+kB.*T.*log(X10D))+X01D.*(E01+kB.*T.*log(X01D)))+ ...,
    (X00E.*(E00+kB.*T.*log(X00E))+X11E.*(E11+kB.*T.*log(X11E))+X10E.*(E10+kB.*T.*log(X10E))+X01E.*(E01+kB.*T.*log(X01E)))+ ...,
    (X00F.*(E00+kB.*T.*log(X00F))+X11F.*(E11+kB.*T.*log(X11F))+X10F.*(E10+kB.*T.*log(X10F))+X01F.*(E01+kB.*T.*log(X01F)));


FF1=(X1.*log(X1)+(1-X1).*log(1-X1)+X2.*log(X2)+(1-X2).*log(1-X2)+X3.*log(X3)+(1-X3).*log(1-X3)+X4.*log(X4)+(1-X4).*log(1-X4))./4;


F=2*N*kB.*T.*(X1.*L1X+X2.*L2X+X3.*L3X+X4.*L4X-log(part4))-N*FF2+5*N*kB.*T.*FF1;
E4=(X1111*E1111+X1110*E1110+X1101*E1101+X1011*E1011+X0111*E0111+X1100*E1100+E1010*X1010+E1001*X1001+X0011*E0011+X0101*E0101+X0110*E0110+X1000*E1000+X0100*E0100+X0010*E0010+X0001*E0001+X0000*E0000)*N;
E2=X00A*E00+X11A*E11+X10A*E10+X01A*E01+X00B*E00+X11B*E11+X10B*E10+X01B*E01+X00C*E00+X11C*E11+X10C*E10+X01C*E01+X00D*E00+X11D*E11+X10D*E10+X01D*E01+X00E*E00+X11E*E11+X10E*E10+X01E*E01+X00F*E00+X11F*E11+X10F*E10+X01F*E01;
EE=(2*E4-N*E2);

%

S=-(F-EE)./T;

RF=F/kB/N;
RE=EE/kB/N;
RS=S/kB/N;

XT=(X1+X2+X3+X4)./4;

    for j=1:1501
        if abs(XT(j,1)-x)<0.0002
        if F(j,1)<minF(j,1)
            minF(j,1)=F(j,1);
            FFF(j,1)=RF(j,1);
            SSS(j,1)=RS(j,1);
            EEE(j,1)=RE(j,1);
            parameter1(j,1)=xxx;
            parameter2(j,1)=XT(j,1);
        end
        end
    end


toc




end
j=0;
for i=1:1500
    C(i)=(EEE(i+1)-EEE(i))/(T(i+1)-T(i));
    if C(i)<10
        j=j+1;
        Cv(j)=C(i);
        TTT(j)=(T(i)+T(i+1))/2;
    end
end
figure(1)
plot(TTT,Cv)
figure(2)
 plot(T,SSS)
figure(3)
 plot(T,EEE)
figure(4)
 plot(T,FFF)