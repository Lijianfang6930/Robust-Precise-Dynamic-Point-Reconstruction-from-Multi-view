function [A,b] = buildAb(uv,M)

% 生成一组3D重建时用的A和b
 %%% 输入：uv,M
 % uv    一行       为输入的多个镜头的2D点数据
 %                  格式：每两列为一个镜头的数据，m个镜头，有2m列  
 % M    camN*11     为输入的多个镜头的M矩阵
 %                  格式：每一行为一个镜头的M矩阵，n个镜头，为n*11的矩阵
 %%% 输出：
 % A    2camN*3     每2行对应一个镜头的uv
 % b    2camN*1     每2行对应一个镜头的uv
 
 
 
 M_l=size(M,1);   %求取M的行数
 uv = uv(:) ;
  
 A=zeros(2*M_l,3); b=zeros(2*M_l,1);%不定长A和b
 for icam = 1:M_l
     A(icam*2-1,:) = uv(icam*2-1) * M(icam,9:11) - M(icam,1:3) ;
     A(icam*2,:)   = uv(icam*2)   * M(icam,9:11) - M(icam,5:7)  ;
     b(icam*2-1)   = M(icam,4)-uv(icam*2-1) ;
     b(icam*2)     = M(icam,8)-uv(icam*2) ;
 end
 
