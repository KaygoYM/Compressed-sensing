function x=SP_paper(Phi,y,K)
%SP�㷨

%��ȡPhi���������������
[M,N]=size(Phi);

%��ʼ������
%��Phi��ÿ����y����أ��õ�һ��N*1�ľ���(������)
correlation=Phi'*y;
%��correlationȡ����ֵ�����򣬰��Ӵ�С��˳��
[var,pos] = sort(abs(correlation),'descend');
%����һ���ռ�T�����ڼ�¼Phi��������ֵ
T=[];T=union(T,pos(1:K));
y_r=resid_paper(y,Phi(:,T));

%����
%ʹ��������ʽ��do---while�ṹ
% while(1)
%	 if(condition)
%		 break;
%	 end
% end
count=1;
while(1)
  %���ݲв��������ӵ��������õ�T_add
  correlation=Phi'*y_r;
  [var,pos] = sort(abs(correlation),'descend');
  T_add=union([],pos(1:K));
  
  %�ϲ����е�T��T_add
  T=union(T,T_add);
  
  %
  x_p=((Phi(:,T)'*Phi(:,T))\eye(length(T)))*Phi(:,T)'*y;%proj_paper(y,Phi(:,T));
  %�����±��¼T
  [var,pos] = sort(abs(x_p),'descend');
  %ȡǰK�����ֵ
  T=union([],T(pos(1:K)));
  
  %�����µĲв�
  y_r_n=resid_paper(y,Phi(:,T));
  
  %�ж��Ƿ��˳�ѭ��,����Ϊ������100��
  if(norm(y_r_n)>=norm(y_r) || count>100)
    break;
  end
  
  %�����˳�ѭ��,������һ�ֵĵ���
  y_r=y_r_n;
  
  count=count+1;
end

%�˳�ѭ�����������������
x=zeros(N,1);
x(T)=((Phi(:,T)'*Phi(:,T))\eye(length(T)))*Phi(:,T)'*y;


end





function y_r=resid_paper(y,Phi)
%����y��Phi�ϵ�ͶӰ�в�

%��ȡ����Phi������������,Mû����
[M,N]=size(Phi);

%�жϾ���(Phi'*Phi)�Ƿ����
if(rank(Phi'*Phi)~=N)
  error('���󲻿���');
end

y_p=Phi*((Phi'*Phi)\eye(N))*Phi'*y;
y_r=y-y_p;
end