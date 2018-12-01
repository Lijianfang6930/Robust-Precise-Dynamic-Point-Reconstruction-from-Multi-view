% ======================================================================
% Step3.ƥ��
% ======================================================================

% ==== Step3.1 �������ü����ݶ�ȡ ====
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
%k = 4; %��ֵѡȡ����������buildF(k),kֵԽ����ֵԽ��
mincam = 3 ; %ƥ���ؽ�����ʹ�þ�ͷ����

% usecam = [ 7 8 9 11 12 13 14 15 ]  %ʹ�õľ�ͷ��ţ�������������������������������������������������
usecam = 1:camN ;
usecamN = length(usecam) ;

% % ��ȡ.vc�ļ�
% filename = '100Marker1' ;
% % datapath = 'input\Nokov_Yanjiao_20160725\VCFiles\' ;
% datapath ='E:\20161018 չ������\NOKOV_Test_All\Nokov_Test_20161018_Wand5\VCFiles\' ;
% usecam = 1:camN ;
% fprintf('\n��ȡ��ת��%s���ݣ�\n',filename)
% frameN = length(usecam) ; %ÿ����ͷ��֡��
% uv = cell(1,length(usecam)) ;
% i=1 ;
% for icam = usecam
% uv{i} = readvc([datapath,filename,'\',filename,'.vc',int2str(icam)]) ;
% frameN(i) =length(uv{i}) ; 
% fprintf('��%d�ž�ͷ���ݶ�ȡת�����\n',icam)
% i=i+1 ;
% end
% uvc = uv ;
load input\uvc_100Marker1

% ��ȡ�����ļ�
datapath = 'input\' ;
M0 = importdata([datapath,'M1.txt']) ;
inexI_hat = importdata([datapath,'inexI.txt']) ;
kk = importdata([datapath,'kk.txt']) ;
uvc = uvd;
% ֻ��ȡ��ǰ����ʹ�õľ�ͷ������
%uvc = readvc(G:\lijianfang\�����Ƽ�\Code\input\20161018 չ������\NOKOV_Test_All\Nokov_Test_20161018_Wand5\VCFiles\50Marker1\50Marker1.vc)
uvc = uvc(usecam) ;
M = M0(usecam,:);
inexI_hat = inexI_hat(usecam,:) ;
kk = kk(usecam,:) ;

%�����������Fundation
F = nan(3*usecamN) ;  %Fundation
for i = 1:usecamN-1
    for j = i+1:usecamN
        [F(3*i-2:3*i,3*j-2:3*j), F(3*j-2:3*j,3*i-2:3*i)] = MtoF(M(i,:),M(j,:)) ; %[F12,F21]
    end
end


% ��֡�Ŵ洢uv�����ݽṹ��������֡��ƥ�����
uvf = []; %uv_frame 
frameN = length(uvc{1}) ;
for icam = 1:length(uvc) %i_cam
    for iframe = 1:length(uvc{icam}) %i_frame
        uvf{iframe}{icam} = uvc{icam}{iframe} ;
    end
end

% ==== Step3.2 ƥ�䡢�ؽ� ====

clear match_frame
xyzn=cell(length(uvf),1); %ÿһ֡�����е��3D����
markerN = []; %ÿ֡ƥ����ĵ���

% Num_dengyu{k} = [];
% Num_xiaoyu{k} = [];
% Num_dayu{k} = [];


for iframe = 1:1:length(uvf)
%     iframe = 161
    uv = uvf{iframe}; %��֡��uv����
    
    % �������
    for iuv = 1:length(uv)
        uv{iuv} = adddistortion(uv{iuv}, inexI_hat(iuv,:), kk(iuv,:)) ;
    end

    % ƥ��
    match_frame(iframe) = match_v11_0_F3(uv,F,k,mincam) ;
    match_frame(iframe).i = iframe ;
    
    % �ؽ�
    xyz = RebuildMatchframe2(match_frame(iframe),M,30) ;
    xyzn{iframe} = xyz ;
    
    
    
    % ��ͼ
    markerN(iframe) = length(xyz)/3 ; %ÿ֡�ĵ���
    
    %ͳ��
   % if 200*k<length(markerN)<=(k+1)*200
%     if markerN(iframe) == 100
%         Num_dengyu(t) = Num_dengyu(t) + 1;
%         %fprintf('marker�����100֡����%d\n',Num_dengyu);
%         else
%         if markerN(iframe) <= 100
%           Num_xiaoyu(t) = Num_xiaoyu(t) + 1;
%           % fprintf('marker��С��100֡����%d\n',Num_xiaoyu);
%         else
%           Num_dayu(t) = Num_dayu(t) + 1;
%           %fprintf('marker�����100֡����%d\n',Num_dayu);
%         end
%     end
    
    
    fprintf('֡�ţ�%d��������%d\n',iframe,markerN(iframe));
    plot3(xyz(3*(1:markerN(iframe))-2),xyz(3*(1:markerN(iframe))-1),xyz(3*(1:markerN(iframe))),'.') 
    axis([-2300 2300 -1800 1800 -1000 2000]); % ������������ָ��������
    n = length(xyz)/3 ;
    a = [int2str(iframe),' (',int2str(n),')'] ;
    title(a) ;
    % grid on;  
    pause(0.00000000001)
    
end
save('xyz1.mat','xyzn')
% save(['input\xyzn_v11_0_',filename,'_k',int2str(k),'c',int2str(mincam)],'xyzn') ;
% disp('xyz�ļ��ѱ���~')
% 
% % ==== Step3.3 ��ͼ ====
% figure
% while 1
% for i = 1:length(uvf)
% %     for i = 960:1000
% xyz = xyzn{i} ;
% if ~isempty(xyz)
% n = length(xyz)/3 ;
% markerN(i) = n; %ÿ֡�ĵ���
% plot3(xyz(3*(1:n)-2),xyz(3*(1:n)-1),xyz(3*(1:n)),'.' ,'markersize',8)
% axis([-2300 2300 -1800 1800 -1000 2000]); % ������������ָ��������
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
%         fprintf('marker�����100֡����%d\n',Num_dengyu);
%     else
%         if markerN(i) <= 100
%           Num_xiaoyu = Num_xiaoyu + 1;
%            fprintf('marker��С��100֡����%d\n',Num_xiaoyu);
%         else
%           Num_dayu = Num_dayu + 1;
%           fprintf('marker�����100֡����%d\n',Num_dayu);
%         end
%     end
%     end 
% end
% 







