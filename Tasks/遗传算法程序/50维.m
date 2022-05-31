function ga
%遗传算法极小值问题
clear
clc
bn=50;     %个体串长度
inn=50;    %初始种群大小
gnmax=600; %最大迭代次数
pc=0.8;   %交叉概率
pm=0.1;   %变异概率
%产生初始种群

fps=10.24*(rand(inn,bn)-0.5*ones(inn,bn));  %十进制编码，使初始解在[-5.12，5.12]范围内
select=fps;
gn=0;
%第一代解
n = 30;
fr = 4000;
while gn<gnmax,
%分别将种群代入目标函数计算适应度
q = 0;
p = 0;
for i = 1:50,
   for j = 1:50,
       x(j)=select(i,j:j);
       q = q+x(j)^2+10;
       p = p+10*cos(2*pi*x(j));
   end
       y(i)=q-p;
       q=0;
       p=0;    
end
[pmin,i]=min(y);%将当前种群中最优解保存在pmin
pmin1=pmin;
vari=mean(select(i,:));%计算最优染色体的平均值并送给vari，变异自适应调节域值
l=i;%将最优染色体行坐标传送给l
[pmax,z]=max(y);%求取适应度最差的染色体
select(z,:)=select(l,:);%将最优染色体引入种群参与竞争
%进行交叉，进一步选择优秀的解
for i=1:2:50,
    pcc=pro(pc);%根据交叉概率判断是否进行交叉
    if pcc==1,
        if 0<=l-i<=1,%判断是不是最优染色体，如果是不进行交叉
     cross(i,:)=select(i,:);
     cross(i+1,:)=select(i+1,:);
        else
    %将30个变量分为五个范围，相临两个染色体分别从五位进行交叉
    crb1=round(rand*(bn/5-1))+1;  %在[1,6]范围内随机产生一个交叉位
    cross(i,crb1:crb1)=0.1*select(i,crb1:crb1)+0.9*select(i+1,crb1:crb1);
    for j=1:9
   cross(i+1,j*crb1:j*crb1)=0.1*select(i,j*crb1:j*crb1)+0.9*select(i+1,j*crb1:j*crb1);
    cross(i,(j+1)*crb1:(j+1)*crb1)=0.1*select(i,(j+1)*crb1:(j+1)*crb1)+0.9*select(i+1,(j+1)*crb1:(j+1)*crb1);
   
    end
      %  %%%%
   %  cross(i+1,crb1:crb1)=0.1*select(i,crb1:crb1)+0.9*select(i+1,crb1:crb1);
     %cross(i,2*crb1:2*crb1)=0.1*select(i,2*crb1:2*crb1)+0.9*select(i+1,2*crb1:2*crb1);
     %cross(i+1,2*crb1:2*crb1)=0.1*select(i,2*crb1:2*crb1)+0.9*select(i+1,2*crb1:2*crb1);
     %cross(i,3*crb1:3*crb1)=0.1*select(i,3*crb1:3*crb1)+0.9*select(i+1,3*crb1:3*crb1);
     %cross(i+1,3*crb1:3*crb1)=0.1*select(i,3*crb1:3*crb1)+0.9*select(i+1,3*crb1:3*crb1);
     %cross(i,4*crb1:4*crb1)=0.1*select(i,4*crb1:4*crb1)+0.9*select(i+1,4*crb1:4*crb1);
     %cross(i+1,4*crb1:4*crb1)=0.1*select(i,4*crb1:4*crb1)+0.9*select(i+1,4*crb1:4*crb1);  
  %  cross(i,5*crb1:5*crb1)=0.1*select(i,5*crb1:5*crb1)+0.9*select(i+1,5*crb1:5*crb1);
   %  cross(i+1,5*crb1:5*crb1)=0.1*select(i,5*crb1:5*crb1)+0.9*select(i+1,5*crb1:5*crb1); 
   %  cross(i,6*crb1:6*crb1)=0.1*select(i,6*crb1:6*crb1)+0.9*select(i+1,6*crb1:6*crb1);
    % cross(i+1,6*crb1:6*crb1)=0.1*select(i,6*crb1:6*crb1)+0.9*select(i+1,6*crb1:6*crb1);  
    % cross(i,7*crb1:7*crb1)=0.1*select(i,7*crb1:7*crb1)+0.9*select(i+1,7*crb1:7*crb1);
    % cross(i+1,7*crb1:7*crb1)=0.1*select(i,7*crb1:7*crb1)+0.9*select(i+1,7*crb1:7*crb1);
    % cross(i,8*crb1:8*crb1)=0.1*select(i,8*crb1:8*crb1)+0.9*select(i+1,8*crb1:8*crb1);
   %  cross(i+1,8*crb1:8*crb1)=0.1*select(i,8*crb1:8*crb1)+0.9*select(i+1,8*crb1:8*crb1);
  %   cross(i,9*crb1:9*crb1)=0.1*select(i,9*crb1:9*crb1)+0.9*select(i+1,9*crb1:9*crb1);
    %cross(i+1,9*crb1:9*crb1)=0.1*select(i,9*crb1:9*crb1)+0.9*select(i+1,9*crb1:9*crb1);


        end
    else
    cross(i,:)=select(i,:);
    cross(i+1,:)=select(i+1,:);
  end
end
%变异操作
for i=1:50;
pmm=pro(pm);
if pmm==1,
    if i==l,
      change(i,:)=cross(i,:);
    else
    change(i,:)=cross(i,:);
    chb=round(rand*(bn-1))+1;  %在[1,bn]范围内随机产生一个变异位
    change(i,chb:chb)= vari*(rand-0.5);%变异为变化为一个[-5.12，5.12]的随机数
     end
else
      change(i,:)=cross(i,:);
end 
end
select=change;
gn=gn+1;
t=x
plot(gn,pmin,'g.')
title('历代适应度变化','fonts',10);
xlabel('1-600代');
ylabel('每代的最优个体');
hold on
end
gn
pmin%最优解
strt=['最大迭代代数' num2str(gn)];
text(390,600,strt);
strt=['最小值' num2str(pmin)];
text(390,500,strt);

%**********************************************************************%



%根据交叉概率判断是否交叉函数
function pcc=pro(pc);
test(1:100)=0;
b=round(100*pc);
test(1:b)=1;
n=round(rand*99)+1;
pcc=test(n);   
end
end

 