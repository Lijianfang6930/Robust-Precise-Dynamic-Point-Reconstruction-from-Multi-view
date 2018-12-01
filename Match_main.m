% ======================================================================
% Step3.匹配
% ======================================================================

% ==== Step3.1 参数设置及数据读取 ====
clear all
close all
%  Num_dengyu =[];
%  Num_xiaoyu =[];
%  Num_dayu = [];
 k =3
%     t = floor((k-3)*10);
% 
%     Num_dengyu(t) = 0;
%     Num_xiaoyu(t) = 0;
%      Num_dayu(t) = 0;
camN = 16 ;
%k = 4; %阈值选取参数，用于buildF(k),k值越大阈值越大
mincam = 3 ; %匹配重建最少使用镜头个数

% usecam = [ 7 8 9 11 12 13 14 15 ]  %使用的镜头编号！！！！！！！！！！！！！！！！！！！！！！！！！
usecam = 1:camN ;
usecamN = length(usecam) ;

% % 读取.vc文件
% filename = '100Marker1' ;
% % datapath = 'input\Nokov_Yanjiao_20160725\VCFiles\' ;
% datapath ='E:\20161018 展会数据\NOKOV_Test_All\Nokov_Test_20161018_Wand5\VCFiles\' ;
% usecam = 1:camN ;
% fprintf('\n读取并转换%s数据：\n',filename)
% frameN = length(usecam) ; %每个镜头的帧数
% uv = cell(1,length(usecam)) ;
% i=1 ;
% for icam = usecam
% uv{i} = readvc([datapath,filename,'\',filename,'.vc',int2str(icam)]) ;
% frameN(i) =length(uv{i}) ; 
% fprintf('第%d号镜头数据读取转换完毕\n',icam)
% i=i+1 ;
% end
% uvc = uv ;
load input\uvc_100Marker1

% 读取参数文件
datapath = 'input\' ;
M0 = importdata([datapath,'M1.txt']) ;
inexI_hat = importdata([datapath,'inexI.txt']) ;
kk = importdata([datapath,'kk.txt']) ;
uvc = uvd;
% 只读取当前正在使用的镜头的数据
%uvc = readvc(G:\lijianfang\度量科技\Code\input\20161018 展会数据\NOKOV_Test_All\Nokov_Test_20161018_Wand5\VCFiles\50Marker1\50Marker1.vc)
uvc = uvc(usecam) ;
M = M0(usecam,:);
inexI_hat = inexI_hat(usecam,:) ;
kk = kk(usecam,:) ;

%计算基础矩阵Fundation
F = nan(3*usecamN) ;  %Fundation
for i = 1:usecamN-1
    for j = i+1:usecamN
        [F(3*i-2:3*i,3*j-2:3*j), F(3*j-2:3*j,3*i-2:3*i)] = MtoF(M(i,:),M(j,:)) ; %[F12,F21]
    end
end


% 按帧号存储uv的数据结构，便于逐帧地匹配操作
uvf = []; %uv_frame 
frameN = length(uvc{1}) ;
for icam = 1:length(uvc) %i_cam
    for iframe = 1:length(uvc{icam}) %i_frame
        uvf{iframe}{icam} = uvc{icam}{iframe} ;
    end
end

% ==== Step3.2 匹配、重建 ====

clear match_frame
xyzn=cell(length(uvf),1); %每一帧的所有点的3D坐标
markerN = []; %每帧匹配出的点数

% Num_dengyu{k} = [];
% Num_xiaoyu{k} = [];
% Num_dayu{k} = [];


for iframe = 1:1:length(uvf)
%     iframe = 161
    uv = uvf{iframe}; %单帧的uv数据
    
    % 畸变矫正
    for iuv = 1:length(uv)
        uv{iuv} = adddistortion(uv{iuv}, inexI_hat(iuv,:), kk(iuv,:)) ;
    end

    % 匹配
    match_frame(iframe) = match_v11_0_F3(uv,F,k,mincam) ;
    match_frame(iframe).i = iframe ;
    
    % 重建
    xyz = RebuildMatchframe2(match_frame(iframe),M,30) ;
    xyzn{iframe} = xyz ;
    
    
    
    % 画图
    markerN(iframe) = length(xyz)/3 ; %每帧的点数
    
    %统计
   % if 200*k<length(markerN)<=(k+1)*200
%     if markerN(iframe) == 100
%         Num_dengyu(t) = Num_dengyu(t) + 1;
%         %fprintf('marker点等于100帧数：%d\n',Num_dengyu);
%         else
%         if markerN(iframe) <= 100
%           Num_xiaoyu(t) = Num_xiaoyu(t) + 1;
%           % fprintf('marker点小于100帧数：%d\n',Num_xiaoyu);
%         else
%           Num_dayu(t) = Num_dayu(t) + 1;
%           %fprintf('marker点大于100帧数：%d\n',Num_dayu);
%         end
%     end
    
    
    fprintf('帧号：%d，点数：%d\n',iframe,markerN(iframe));
    plot3(xyz(3*(1:markerN(iframe))-2),xyz(3*(1:markerN(iframe))-1),xyz(3*(1:markerN(iframe))),'.') 
    axis([-2300 2300 -1800 1800 -1000 2000]); % 设置坐标轴在指定的区间
    n = length(xyz)/3 ;
    a = [int2str(iframe),' (',int2str(n),')'] ;
    title(a) ;
    % grid on;  
    pause(0.00000000001)
    
end
save('xyz1.mat','xyzn')
% save(['input\xyzn_v11_0_',filename,'_k',int2str(k),'c',int2str(mincam)],'xyzn') ;
% disp('xyz文件已保存~')
% 
% % ==== Step3.3 画图 ====
% figure
% while 1
% for i = 1:length(uvf)
% %     for i = 960:1000
% xyz = xyzn{i} ;
% if ~isempty(xyz)
% n = length(xyz)/3 ;
% markerN(i) = n; %每帧的点数
% plot3(xyz(3*(1:n)-2),xyz(3*(1:n)-1),xyz(3*(1:n)),'.' ,'markersize',8)
% axis([-2300 2300 -1800 1800 -1000 2000]); % 设置坐标轴在指定的区间
% % axis([-2000 2000 -2000 2000 0 1000]); 
% a = [int2str(i),' (',int2str(n),')'] ;
% title(a) ;
% pause(0.015)
% % pause(0.1)
%     end
% end
% 
% end
% Num_dengyu = 0;
% Num_xiaoyu = 0;
% Num_dayu = 0;
% while length(markerN)<=200
%     for i = 1:200;
%     if markerN(i) == 100
%         Num_dengyu = Num_dengyu + 1;
%         fprintf('marker点等于100帧数：%d\n',Num_dengyu);
%     else
%         if markerN(i) <= 100
%           Num_xiaoyu = Num_xiaoyu + 1;
%            fprintf('marker点小于100帧数：%d\n',Num_xiaoyu);
%         else
%           Num_dayu = Num_dayu + 1;
%           fprintf('marker点大于100帧数：%d\n',Num_dayu);
%         end
%     end
%     end 
% end
% 







