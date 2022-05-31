%用遗传算法进行简单函数的优化,显示中间过程
function ga
clear

bn=24; %个体串长度1.024*10^7<2^24
inn=50; %初始种群大小
gnmax=200;  %最大代数
pc=0.8; %交叉概率
pm=0.05; %变异概率

%产生初始种群
s=round(rand(inn,bn));

gnf1=5;
gnf2=20;

%计算适应度,返回适应度f和累积概率p
[f,p]=objf(s);  

gn=1;
while gn<gnmax+1
    xp=-5.12:0.01:5.12;
    yp=ft(xp);
    for d=1:inn
        xi=n2to10(s(d,:));
        xdi(d)=-5.12+xi*10.24/(power(2,bn)-1);
    end
    yi=ft(xdi);
    plot(xp,yp,'b-',xdi,yi,'g*');
    strt=['当前代数 gn=' num2str(gn)];
    text(-0.75,1,strt);
    text(-0.75,3.5,'*  当前种群','Color','g');
    if gn<gnf1
       
    end
    hold on;
           
    for j=1:2:inn
      %选择操作
      seln=sel(s,p);
      xs1=n2to10(s(seln(1),:));
      xds1=-5.12+xs1*10.24/(power(2,bn)-1);
      ys1=ft(xds1);
      xs2=n2to10(s(seln(2),:));
      xds2=-5.12+xs2*10.24/(power(2,bn)-1);
      ys2=ft(xds2);
      hold on;
      drawnow;
      plot(xds1,ys1,'r*',xds2,ys2,'r*');
      %交叉操作
      scro=cro(s,seln,pc);
      scnew(j,:)=scro(1,:);
      scnew(j+1,:)=scro(2,:);
      
      %变异操作
      smnew(j,:)=mut(scnew(j,:),pm);
      smnew(j+1,:)=mut(scnew(j+1,:),pm);
      
  end
  drawnow;
  text(-0.7,3.3,'*  选择后','Color','r');
  if gn<gnf1
     
  end
  
  for d=1:inn
      xc=n2to10(scnew(d,:));
      xdc(d)=-5.12+xc*10.24/(power(2,bn)-1);
  end
  yc=ft(xdc);
  drawnow;
  plot(xdc,yc,'m*');
  text(-0.65,3.1,'*  交叉后','Color','m');
  if gn<gnf1
  
  end
  hold on;
  
  for d=1:inn
      xm=n2to10(smnew(d,:));
      xdm(d)=-5.12+xm*10.24/(power(2,bn)-1);
  end
  ym=ft(xdm);
  drawnow;
  plot(xdm,ym,'c*');
  text(-0.60,2.9,'*  变异后','Color','c');
  
  if gn<gnf2
    
  end
  hold off;
  s=smnew;  %产生了新的种群
   
   %计算新种群的适应度   
   [f,p]=objf(s);
   
   %记录当前代最好和平均的适应度
   [fmax,nmax]=max(f);
   fmean=mean(f);
   ymax(gn)=fmax;
   ymean(gn)=fmean;
   
   %记录当前代的最佳个体
   x=n2to10(s(nmax,:));
   xx=-5.12+x*10.24/(power(2,bn)-1);
   xmax(gn)=xx;
   
   gn=gn+1;
end
gn=gn-1;

figure(2);
subplot(2,1,1);
plot(1:gn,[ymax;ymean]);
title('历代适应度变化','fonts',10);
legend('最大适应度','平均适应度');
string1=['最终适应度',num2str(ymax(gn))];
gtext(string1);
subplot(2,1,2);
plot(1:gn,xmax,'r-');
legend('自变量');
string2=['最终自变量',num2str(xmax(gn))];
gtext(string2);
end

%“交叉”操作
function scro=cro(s,seln,pc);

[inn bn]=size(s);

if rand<pc
   chb=ceil(rand*(bn-1));  %在[1,bn-1]范围内随机产生一个交叉位
   scro(1,:)=[s(seln(1),1:chb) s(seln(2),chb+1:bn)];
   scro(2,:)=[s(seln(2),1:chb) s(seln(1),chb+1:bn)];
else
   scro(1,:)=s(seln(1),:);
   scro(2,:)=s(seln(2),:);
end  
end

%目标函数
function y=ft(x);

    y=x.*x-10*cos(2*pi*x)+10;
 
end

%“变异”操作
function snnew=mut(snew,pm);

bn=size(snew,2);
snnew=snew;

if rand<pm
   chb=ceil(rand*bn);  %在[1,bn]范围内随机产生一个变异位
   snnew(chb)=abs(snew(chb)-1);
end   
end

%将2进制数转换为10进制数
function x=n2to10(s);

bn=size(s,2);
x=s(bn);
for i=1:bn-1
   x=x+s(bn-i)*power(2,i);
end
end

%计算适应度函数

function [f,p]=objf(s);

[inn bn]=size(s);  %读取种群大小, 有inn个个体, 个体长度为bn

for i=1:inn
   x=n2to10(s(i,:));  %讲二进制转换为十进制
   xx=-5.12+x*10.24/(power(2,bn)-1);  %转化为[-5.12,5.12]区间的实数
   f(i)=1/(1+ft(xx));  %计算函数值，即适应度
end
f=f';

%计算选择概率
fsum=0;
for i=1:inn
   fsum=fsum+f(i)*f(i);
end
for i=1:inn
   ps(i)=f(i)*f(i)/fsum;
end

%计算累积概率
p(1)=ps(1);
for i=2:inn
   p(i)=p(i-1)+ps(i);
end
p=p';
end

%“选择”操作
function seln=sel(s,p);

inn=size(p,1);

%从种群中选择两个个体
for i=1:2
   r=rand;  %产生一个随机数
   prand=p-r;
   j=1;
   while prand(j)<0
       j=j+1;
   end
   seln(i)=j; %选中个体的序号
end
end
