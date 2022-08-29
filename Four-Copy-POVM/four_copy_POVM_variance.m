%Aug 2022: file for reading the four-copy POVM

POVMS4=csvread('4outcome_POVM.csv');
eps4=csvread('4outcome_est.csv');
%eps4=real(eps4);
 POVMS4*POVMS4';

%decoherence, epsilon in paper
p=0.5;
%angle
th=0
%initial state
psi=[1,0];
%rotated state
psith=Ry(th)*Rx(th)*psi';

rhoth=psith*psith';
%decoehred state
rho=(1-p)*rhoth+p*eye(2)/2;


rho2=kron(rho,rho);
rho3=kron(rho2,rho);
%four-copy state
rho4=kron(rho3,rho);

%variances 
vx=0;
vy=0;
for hh=1:16
    Pi4{hh}=POVMS4(:,hh)*POVMS4(:,hh)';
    prob4(hh)=trace(rho4*Pi4{hh});

    vx=vx+prob4(hh)*eps4(1,hh)^2;
    vy=vy+prob4(hh)*eps4(2,hh)^2;
end

%comparison of this variance with the variances in the paper
variance4copy=4*(vx+vy)
hol=(4-2*p)/((1-p)^2)
nag1=4/(1-p)^2
nag2=2*(4-2*p+p^2)/(2*(1-p)^2)
nag3=12.716


