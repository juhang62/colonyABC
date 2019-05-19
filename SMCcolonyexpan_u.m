clear
load('elife colony data processed.mat')
m = 1; %cylidrical
x = linspace(0,700,701);
t = [16.2000   20.7000   24.6000   28.5000];
%t = tobs - tobs(1);
%t = [0 1 2 3];
resume_flag=0;
sig_q=1;
if resume_flag==1
    load('resume.mat')
    epss=ep*0.9.^(1:20);
    wnew=w;
else
    a=1;b=1;D=1;
    theta=[a;b;D];
    epss=3e6*0.7.^(1:150);
    npar=500;
    thetapop=zeros(3,npar);
    w=ones(npar,1);
    wnew=w;
    thetapopnew=thetapop;
end

for ieps=1:length(epss)
    %ieps
    ep=epss(ieps)
    xini=xo{1};
    yini=yo{1};
    parfor ipar=1:npar
        d=1e10;theta=zeros(3,1);
        rej=0;
        while d>ep || any(theta<=0) || any(theta>=20)
            if ieps==1 && ~resume_flag
                theta=rand(3,1)*20;
            else
                ind=randsample(npar,1,true,w);
                theta=thetapop(:,ind);
                theta=theta+normrnd(0,sig_q,3,1);
                if any(theta<=0) || any(theta>=20)
                    rej=rej+1;
                    continue
                end
            end
            sol = pdepe(m,@(x,t,u,DuDx) pdefun(x,t,u,DuDx,theta),@(x) pdeic(x,xini,yini),@pdebc,x,t);
            d = dist(x,sol,xo,yo);
            if d>ep
                rej=rej+1;
                if rej>200
                    disp('fuck')
		    quit
                end
                continue
            end
            thetapopnew(:,ipar)=theta;
            if ieps>1 || resume_flag 
                Ktt=exp(-sum((thetapop-repmat(theta,[1 npar])).^2));
                wnew(ipar)=1/sum(Ktt*w);
            end
        end
%         rej
%         ipar
    end
    thetapop=thetapopnew;
    w=wnew/(sum(wnew));
    save('resume','thetapop','w','ep','npar')
end



%%
% for i=1:length(t)
%     plot(x,sol(i,:))
%     title(['time=', num2str(t(i))]);
%     xlabel('radius')
%     ylabel('height')
%     pause(0.1)
% end


function [c,f,s] = pdefun(x,t,u,DuDx,theta)
c = 1;
a=theta(1);
b=theta(2);
D=theta(3);
uc=200;
f =uc * D* (u/uc)^2 .* DuDx;
%s = [u(2)*tanh(u(1))-0.1*u(1); -u(2)*tanh(u(1)); 0.1*u(1)];
s = b*tanh(a*u);
end

function u0 = pdeic(x,xobs,yobs)
%u0 = [0.1*exp(-x^2);1];
if x<max(xobs) && x>min(xobs)
    u0=interp1(xobs,yobs,x);
else
    u0=0;
end
end

function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t)
pl = [0];
ql = [1];
pr = [0];
qr = [1];
end

function out=dist(x,sol,xobs,yobs)
nobs=length(xobs);
out=0;
for i=2:nobs
    ypred=interp1(x,sol(i,:),xobs{i});
    out=out+sum((ypred-yobs{i}).^2);
end
end
