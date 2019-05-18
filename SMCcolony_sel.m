clear
load('elife colony data processed.mat')
m = 1; %cylidrical
x = linspace(0,800,801);
t = [16.2000   20.7000   24.6000   28.5000];
%t = tobs - tobs(1);
%t = [0 1 2 3];
resume_flag=1;
sig_q=1;
if resume_flag==1
    load('resume.mat')
    epss=ep*0.85.^(1:3);
    wnew=w;
    mdlpopnew=mdlpop;
else
    a=1;b=1;D=1;
    theta=[a;b;D];
    %epss=3e6*0.7.^(1:150);
    epss=linspace(3e6,1e5,8);
    npar=1024;
    thetapop=zeros(3,npar);
    w=ones(npar,1);
    wnew=w;
    thetapopnew=thetapop;
    mdlpop=zeros(npar,1);
    mdlpopnew=mdlpop;
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
            mdl=randsample(2,1);
            wmdl=w.*double(mdlpop==mdl);
            if ieps==1 && ~resume_flag
                theta=rand(3,1)*20;
            else
                ind=randsample(npar,1,true,wmdl);
                theta=thetapop(:,ind);
                theta=theta+normrnd(0,sig_q,3,1);
                if any(theta<=0) || any(theta>=20)
                    rej=rej+1;
                    continue
                end
            end
            sol = pdepe(m,@(x,t,u,DuDx) pdefun(x,t,u,DuDx,theta,mdl),@(x) pdeic(x,xini,yini),@pdebc,x,t);
            d = dist(x,sol,xo,yo);
            if d>ep
                rej=rej+1;
                if rej>300
                    error('fuck')
                end
                continue
            end
            mdlpopnew(ipar)=mdl;
            thetapopnew(:,ipar)=theta;
            if ieps>1 || resume_flag 
                Ktt=exp(-sum((thetapop-repmat(theta,[1 npar])).^2));
                wnew(ipar)=1/sum(Ktt*wmdl); %w or wmdl
            end
        end
%         rej
%         ipar
    end
    thetapop=thetapopnew;
    mdlpop=mdlpopnew;
    ind1 = (mdlpop==1);
    wnew(ind1) = wnew(ind1)/(sum(wnew(ind1)));
    wnew(~ind1)= wnew(~ind1)/(sum(wnew(~ind1)));
    w=wnew;
    fname=strcat('resume-ep',num2str(round(epss(ieps)),'%d'));
    save(fname,'thetapop','mdlpop','w','ep','npar')
end



%%
% for i=1:length(t)
%     plot(x,sol(i,:))
%     title(['time=', num2str(t(i))]);
%     xlabel('radius')
%     ylabel('height')
%     pause(0.1)
% end


function [c,f,s] = pdefun(x,t,u,DuDx,theta,mdl)
c = 1;
a=theta(1);
b=theta(2);
D=theta(3);
if mdl==1
    f = D*DuDx;
elseif mdl==2
    f=D*u*DuDx;
end
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
