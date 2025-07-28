absVel = [0,10,20,40,80,160,320,640];
relVel = [0,1/4,1/2,1,2,4,NaN,NaN];

respVel40Abs  = NaN(6,size(absVel,2));
respVel80Abs  = NaN(6,size(absVel,2));
respVel160Abs = NaN(6,size(absVel,2));
respVel40Rel  = NaN(6,size(absVel,2));
respVel80Rel  = NaN(6,size(absVel,2));
respVel160Rel = NaN(6,size(absVel,2));

resp0 = 3;
resp0_25 = 2.6;
resp0_5 = 2.3;
resp1 = 2;
resp2 = 0.5;
resp4 = 1.8;

respVal40Gen  = [resp0,resp0_25,resp0_5,resp1,resp2,resp4,NaN,NaN];
respVal80Gen  = [resp0,NaN,resp0_25,resp0_5,resp1,resp2,resp4,NaN];
respVal160Gen = [resp0,NaN,NaN,resp0_25,resp0_5,resp1,resp2,resp4];
respValRel = [resp0,resp0_25,resp0_5,resp1,resp2,resp4,NaN,NaN];

for i = 1:size(absVel,2)
    respVel40Abs(:,i)  = respVal40Gen(i)  + 0.5*(-1).^(round(rand(1,6))).*(rand(1,6));
    respVel80Abs(:,i)  = respVal80Gen(i)  + 0.5*(-1).^(round(rand(1,6))).*(rand(1,6));
    respVel160Abs(:,i) = respVal160Gen(i) + 0.5*(-1).^(round(rand(1,6))).*(rand(1,6));
    respVel40Rel(:,i)  = respValRel(i) + 0.5*(-1).^(round(rand(1,6))).*(rand(1,6));
    respVel80Rel(:,i)  = respValRel(i) + 0.5*(-1).^(round(rand(1,6))).*(rand(1,6));
    respVel160Rel(:,i) = respValRel(i) + 0.5*(-1).^(round(rand(1,6))).*(rand(1,6));
end

x_rel = reshape(repmat(relVel,[6,3]),[],1);
y_rel = reshape([respVel40Rel,respVel80Rel,respVel160Rel],[],1);

x_abs = reshape(repmat(absVel,[6,3]),[],1);
y_abs = reshape([respVel40Abs,respVel80Abs,respVel160Abs],[],1);

[b_rel,dev_rel,stats_rel] = glmfit(x_rel,y_rel,'gamma');
[b_abs,dev_abs,stats_abs] = glmfit(x_abs,y_abs,'gamma');

yfit_rel = glmval(b_rel,x_rel,'reciprocal');
yfit_abs = glmval(b_abs,x_abs,'reciprocal');

figure
plot(x_rel,y_rel,'o',x_rel,yfit_rel,'-');

figure
plot(x_abs,y_abs,'o',x_abs,yfit_abs,'-');