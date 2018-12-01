function [F12,F21]=MtoF(M1,M2)

% 函数功能：输入M1，M2，输出基础矩阵F12，F21% % % % % % % % % % % % % % % % 
%其中M1，M2，是1*11 行的格式

% % % % % 格式的转换% % % % 

M=zeros(11,2);
M(:,1)=M1';
M(:,2)=M2';



%读取txt文件begin
% M=importdata('2M.txt');%读取两个镜头的M值,第一列为1号镜头，第二列为2号镜头
%********镜头1的M31、m41********
M31(1,1:3)=M(1:3,1); 
M31(2,1:3)=M(5:7,1);
M31(3,1:3)=M(9:11,1);
m41(1:2,1)=M(4:4:8,1);
m41(3,1)=1;
%********镜头2的M32、m42********
M32(1,1:3)=M(1:3,2); 
M32(2,1:3)=M(5:7,2);
M32(3,1:3)=M(9:11,2);
m42(1:2,1)=M(4:4:8,2);
m42(3,1)=1;

%********计算F21********
% m=m42-M32*inv(M31)*m41;
% mx=[0   -m(3,1)   m(2,1)
%     m(3,1)   0   -m(1,1)
%     -m(2,1)   m(1,1)   0];
% 
% F21=mx*M32*inv(M31);
% F21=F21/F21(3,3);
% % dlmwrite('F21_out.txt',F21,'delimiter','\t', 'precision', 10);

%********计算F12********

% m=m41-M31*inv(M32)*m42;
m=m41-M31*(M32\m42);
mx=[0   -m(3,1)   m(2,1)
    m(3,1)   0   -m(1,1)
    -m(2,1)   m(1,1)   0];

% F12=mx*M31*inv(M32);
F12=mx*M31/M32;
F12=F12/F12(3,3);

% dlmwrite('F12_out.txt',F12,'delimiter','\t', 'precision', 10);

F21 = F12' ;

%检验F12与F21的值
% I1=[uv1(1,1:2),1]
% I2=[uv2(1,1:2),1]
% xx12=I1*F12*I2'
% xx21=I2*F21*I1'

end