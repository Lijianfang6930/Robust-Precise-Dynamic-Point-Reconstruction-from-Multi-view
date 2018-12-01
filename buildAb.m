function [A,b] = buildAb(uv,M)

% ����һ��3D�ؽ�ʱ�õ�A��b
 %%% ���룺uv,M
 % uv    һ��       Ϊ����Ķ����ͷ��2D������
 %                  ��ʽ��ÿ����Ϊһ����ͷ�����ݣ�m����ͷ����2m��  
 % M    camN*11     Ϊ����Ķ����ͷ��M����
 %                  ��ʽ��ÿһ��Ϊһ����ͷ��M����n����ͷ��Ϊn*11�ľ���
 %%% �����
 % A    2camN*3     ÿ2�ж�Ӧһ����ͷ��uv
 % b    2camN*1     ÿ2�ж�Ӧһ����ͷ��uv
 
 
 
 M_l=size(M,1);   %��ȡM������
 uv = uv(:) ;
  
 A=zeros(2*M_l,3); b=zeros(2*M_l,1);%������A��b
 for icam = 1:M_l
     A(icam*2-1,:) = uv(icam*2-1) * M(icam,9:11) - M(icam,1:3) ;
     A(icam*2,:)   = uv(icam*2)   * M(icam,9:11) - M(icam,5:7)  ;
     b(icam*2-1)   = M(icam,4)-uv(icam*2-1) ;
     b(icam*2)     = M(icam,8)-uv(icam*2) ;
 end
 
