function ga
%�Ŵ��㷨��Сֵ����
clear
clc
bn=50;     %���崮����
inn=50;    %��ʼ��Ⱥ��С
gnmax=600; %����������
pc=0.8;   %�������
pm=0.1;   %�������
%������ʼ��Ⱥ

fps=10.24*(rand(inn,bn)-0.5*ones(inn,bn));  %ʮ���Ʊ��룬ʹ��ʼ����[-5.12��5.12]��Χ��
select=fps;
gn=0;
%��һ����
n = 30;
fr = 4000;
while gn<gnmax,
%�ֱ���Ⱥ����Ŀ�꺯��������Ӧ��
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
[pmin,i]=min(y);%����ǰ��Ⱥ�����ŽⱣ����pmin
pmin1=pmin;
vari=mean(select(i,:));%��������Ⱦɫ���ƽ��ֵ���͸�vari����������Ӧ������ֵ
l=i;%������Ⱦɫ�������괫�͸�l
[pmax,z]=max(y);%��ȡ��Ӧ������Ⱦɫ��
select(z,:)=select(l,:);%������Ⱦɫ��������Ⱥ���뾺��
%���н��棬��һ��ѡ������Ľ�
for i=1:2:50,
    pcc=pro(pc);%���ݽ�������ж��Ƿ���н���
    if pcc==1,
        if 0<=l-i<=1,%�ж��ǲ�������Ⱦɫ�壬����ǲ����н���
     cross(i,:)=select(i,:);
     cross(i+1,:)=select(i+1,:);
        else
    %��30��������Ϊ�����Χ����������Ⱦɫ��ֱ����λ���н���
    crb1=round(rand*(bn/5-1))+1;  %��[1,6]��Χ���������һ������λ
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
%�������
for i=1:50;
pmm=pro(pm);
if pmm==1,
    if i==l,
      change(i,:)=cross(i,:);
    else
    change(i,:)=cross(i,:);
    chb=round(rand*(bn-1))+1;  %��[1,bn]��Χ���������һ������λ
    change(i,chb:chb)= vari*(rand-0.5);%����Ϊ�仯Ϊһ��[-5.12��5.12]�������
     end
else
      change(i,:)=cross(i,:);
end 
end
select=change;
gn=gn+1;
t=x
plot(gn,pmin,'g.')
title('������Ӧ�ȱ仯','fonts',10);
xlabel('1-600��');
ylabel('ÿ�������Ÿ���');
hold on
end
gn
pmin%���Ž�
strt=['����������' num2str(gn)];
text(390,600,strt);
strt=['��Сֵ' num2str(pmin)];
text(390,500,strt);

%**********************************************************************%



%���ݽ�������ж��Ƿ񽻲溯��
function pcc=pro(pc);
test(1:100)=0;
b=round(100*pc);
test(1:b)=1;
n=round(rand*99)+1;
pcc=test(n);   
end
end

 