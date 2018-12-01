function xyz = RebuildMatchframe2(match_frame,M,threshold3D)
% 重建
% 两两镜头重建拟合
% v9_4中的重建方法
% 
% match_frame   一帧的数据
% xyz           重建出的Marker点坐标，排成一行，每3个元素为一个点的坐标

np = max(match_frame.pass) ; %匹配出来的Marker点数

% 一个一个匹配组地抛点

k=1 ;
for in = 1:np
    ind = match_frame.pass ;
    ind = find(ind==in) ;  %匹配组中的点的位置
    if isempty(ind), continue ; end
    temp = cptxyz(ind, match_frame, M, threshold3D) ;
    if ~isempty(temp)
        xyz(k*3-2:k*3) = temp ;
        k=k+1 ;
    end
end


end %RebuildMatchframe

function xyz = cptxyz(ind, match_frame,M,threshold3D)

% step0.两两重建计算3D坐标

camline = match_frame.camline ;
uv = match_frame.uv ;
n = length(ind) ;  %匹配组中的点数
group0 = camline(:,ind) ;

uvg = uv(:,ind) ;
uvg = uvg(:) ;
Mg = M(group0(1,:),:) ;
[A,b] = buildAb(uvg,Mg) ;
    
answer = zeros(n*(n-1)/2, 5) ; %两两镜头重建的3D坐标
k=1 ;
for i = 1:n-1
    for j = i+1:n
        answer(k,:) = [i,j,(A([i*2-1:i*2,j*2-1:j*2],:)\b([i*2-1:i*2,j*2-1:j*2]))'] ; %rebulid_3D_UnspecCam
        k = k +1 ;
        % answer( (i-1)*(2*n-i)/2 + (j-i),: ) = [i,j, (A([i*2-1:i*2,j*2-1:j*2],:)\b([i*2-1:i*2,j*2-1:j*2]))' ] ; %rebulid_3D_UnspecCam
    end
end

dis = pdist(answer(:,3:5),'cityblock') ; % distance，3D坐标间的曼哈顿距离

if max(dis)<threshold3D, xyz = mean(answer(:,3:5),1) ; return; end


% step1.寻找距离最小值的两点对应的镜头
[mindis, indmindis] = min(dis) ;
if mindis > threshold3D, xyz = []; return; end %该组中各点之间重建差距较大，整组丢掉。

n_answer = size(answer,1) ;
i_mindis = n_answer - floor( ( sqrt(4*((n_answer-1)^2+(n_answer-1)-2*indmindis)+1)+1 ) /2 ) ;
j_mindis = n_answer + indmindis - ((n_answer-1)+(n_answer-i_mindis))*((n_answer-1)-(n_answer-i_mindis)+1)/2 ;

standpoint = [answer(i_mindis,1:2), answer(j_mindis,1:2) ] ;
uniqstandpoint = unique(standpoint) ;

temp = histc(standpoint,uniqstandpoint) ;
igs = uniqstandpoint(temp==max(temp)) ; %index_group_standard
igs = igs(1) ;
xyz_standard = answer(i_mindis,3:5) ; %基准3D坐标


% step2.组中其它镜头上in 的点，挨个与基准点进行重建，重建的3D点与基准3D点距离小于阈值的镜头中的点纳入匹配组中
xyz_new = [ answer(answer(:,1)==igs,:); answer(answer(:,2)==igs,[2,1,3,4,5] ) ] ;
txyz = bsxfun(@minus, xyz_new(:,3:5), xyz_standard) ; % txyz = xyz_new(:,3:5) - repmat(xyz,size(xyz_new,1),1) ;
idxpass = sum(abs(txyz),2) < threshold3D ;
if sum(idxpass)<1 
    xyz=[]; 
    return; 
else
    icam = [igs; xyz_new( idxpass,2 )] ; 
    if length(icam) < 3
        xyz = [];
    else
    xyz = A([icam*2-1;icam*2],:)\b([icam*2-1;icam*2],:) ;
    end
end

end 


