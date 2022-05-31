%���Ŵ��㷨���м򵥺������Ż�,��ʾ�м����
function ga
clear

bn=24; %���崮����1.024*10^7<2^24
inn=50; %��ʼ��Ⱥ��С
gnmax=200;  %������
pc=0.8; %�������
pm=0.05; %�������

%������ʼ��Ⱥ
s=round(rand(inn,bn));

gnf1=5;
gnf2=20;

%������Ӧ��,������Ӧ��f���ۻ�����p
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
    strt=['��ǰ���� gn=' num2str(gn)];
    text(-0.75,1,strt);
    text(-0.75,3.5,'*  ��ǰ��Ⱥ','Color','g');
    if gn<gnf1
       
    end
    hold on;
           
    for j=1:2:inn
      %ѡ�����
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
      %�������
      scro=cro(s,seln,pc);
      scnew(j,:)=scro(1,:);
      scnew(j+1,:)=scro(2,:);
      
      %�������
      smnew(j,:)=mut(scnew(j,:),pm);
      smnew(j+1,:)=mut(scnew(j+1,:),pm);
      
  end
  drawnow;
  text(-0.7,3.3,'*  ѡ���','Color','r');
  if gn<gnf1
     
  end
  
  for d=1:inn
      xc=n2to10(scnew(d,:));
      xdc(d)=-5.12+xc*10.24/(power(2,bn)-1);
  end
  yc=ft(xdc);
  drawnow;
  plot(xdc,yc,'m*');
  text(-0.65,3.1,'*  �����','Color','m');
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
  text(-0.60,2.9,'*  �����','Color','c');
  
  if gn<gnf2
    
  end
  hold off;
  s=smnew;  %�������µ���Ⱥ
   
   %��������Ⱥ����Ӧ��   
   [f,p]=objf(s);
   
   %��¼��ǰ����ú�ƽ������Ӧ��
   [fmax,nmax]=max(f);
   fmean=mean(f);
   ymax(gn)=fmax;
   ymean(gn)=fmean;
   
   %��¼��ǰ������Ѹ���
   x=n2to10(s(nmax,:));
   xx=-5.12+x*10.24/(power(2,bn)-1);
   xmax(gn)=xx;
   
   gn=gn+1;
end
gn=gn-1;

figure(2);
subplot(2,1,1);
plot(1:gn,[ymax;ymean]);
title('������Ӧ�ȱ仯','fonts',10);
legend('�����Ӧ��','ƽ����Ӧ��');
string1=['������Ӧ��',num2str(ymax(gn))];
gtext(string1);
subplot(2,1,2);
plot(1:gn,xmax,'r-');
legend('�Ա���');
string2=['�����Ա���',num2str(xmax(gn))];
gtext(string2);
end

%�����桱����
function scro=cro(s,seln,pc);

[inn bn]=size(s);

if rand<pc
   chb=ceil(rand*(bn-1));  %��[1,bn-1]��Χ���������һ������λ
   scro(1,:)=[s(seln(1),1:chb) s(seln(2),chb+1:bn)];
   scro(2,:)=[s(seln(2),1:chb) s(seln(1),chb+1:bn)];
else
   scro(1,:)=s(seln(1),:);
   scro(2,:)=s(seln(2),:);
end  
end

%Ŀ�꺯��
function y=ft(x);

    y=x.*x-10*cos(2*pi*x)+10;
 
end

%�����족����
function snnew=mut(snew,pm);

bn=size(snew,2);
snnew=snew;

if rand<pm
   chb=ceil(rand*bn);  %��[1,bn]��Χ���������һ������λ
   snnew(chb)=abs(snew(chb)-1);
end   
end

%��2������ת��Ϊ10������
function x=n2to10(s);

bn=size(s,2);
x=s(bn);
for i=1:bn-1
   x=x+s(bn-i)*power(2,i);
end
end

%������Ӧ�Ⱥ���

function [f,p]=objf(s);

[inn bn]=size(s);  %��ȡ��Ⱥ��С, ��inn������, ���峤��Ϊbn

for i=1:inn
   x=n2to10(s(i,:));  %��������ת��Ϊʮ����
   xx=-5.12+x*10.24/(power(2,bn)-1);  %ת��Ϊ[-5.12,5.12]�����ʵ��
   f(i)=1/(1+ft(xx));  %���㺯��ֵ������Ӧ��
end
f=f';

%����ѡ�����
fsum=0;
for i=1:inn
   fsum=fsum+f(i)*f(i);
end
for i=1:inn
   ps(i)=f(i)*f(i)/fsum;
end

%�����ۻ�����
p(1)=ps(1);
for i=2:inn
   p(i)=p(i-1)+ps(i);
end
p=p';
end

%��ѡ�񡱲���
function seln=sel(s,p);

inn=size(p,1);

%����Ⱥ��ѡ����������
for i=1:2
   r=rand;  %����һ�������
   prand=p-r;
   j=1;
   while prand(j)<0
       j=j+1;
   end
   seln(i)=j; %ѡ�и�������
end
end
