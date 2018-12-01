function match_pass = match_v11_0_F3(uv,F,threshold,mincam)
% 按给定的阈值进行匹配
% 流程为：
% 1.粗匹配（双极线）
% 2.投票法筛选（双极线）
% 3.点池删点
% 4.抛点（点残差）
%
% 输入：
% uv             1*camN    cell型，单帧的所有镜头uv信息（camN为镜头数量，下同）
% M              camN*11   每个镜头的向量型M矩阵
% threshold      1*1       点到极线距离的阈值
% mincam         单个数值   匹配中最少使用的镜头个数
% inexI          camN*10   镜头的内外参数
% 输出：
% match_pass 1*n      cell型，每个元素表示匹配出的一个点
%                     cell中的格式，2*N，第一行为镜头号，第二行为镜头中的点号

% 历史版本：
% match_v1_1    在最后加入了去组内错点的步骤，比较组内重建的差异，去掉重建效果不好的uv点，（好像就是点残差……）
% match_v1_2    抛点使用线残差
% match_v1_3    沿用match_v1_1的点残差，将点残差check_points的step2中，
%               两两重建寻找基准uv点时的寻找顺序改为倒向（感觉并没有什么用……）.
% match_v1_4    重写了check_points，修改了点残差的比较方式.
% match_v1_5    在match_v1_4的基础上，修改整个匹配的流程为：
%               1.粗匹配（双极线）
%               2.组内抛点（点残差）
%               3.检漏（点残差）
%               4.点池删点
%               重建结果，闪点变少了，但是点抖得厉害
%               抖点的原因是阈值选大了，调小阈值，抖得就没那么厉害了
% match_v5      在match_v1_5的基础上
%               将threshold_pass的数据结构变成大矩阵的形式.
%               行和列的元素均是： C1P1 C1P2 ... C1Pn_1 C2P1 C2P2 ... C2Pn_2 C3P1 ... CnP1 CnP2 ... CnPn_n
%               其中CiPj表示第i号镜头的第j个点，n_k表示第k个镜头中点的总数为n_k.
%               threshold_pass的第C1P1行、C3P1列表示第1号镜头的1号点与第3号镜头的1号点之间的双极线方法计算值z.
%               子函数相应地做了改变，对同时匹配到同一个镜头中两个点的情况没有做处理。
% match_v5_1_F2 在match_v5的基础上
%               筛选部分screen中，改为判断点到极线的距离阈值，统一阈值为固定值，不再每两两镜头计算不同的阈值
%               函数名后缀带F2的，均为点到极线的距离阈值
% match_v5_2_F2 在match_v5_1_F2的基础上
%               双极线方法粗匹配的步骤中，应保证基准点pb所在的镜头中，没有其他与pa匹配的点
% match_v5_4_F2 在match_v5_2_F2的基础上
%               match_v5_2_F2中第3步检漏的开始部分，会删除已匹配的点；
%               在match_v5_4_F2中，将删除点的步骤放到第2步的开始部分，防止第3步中重新找回.
% match_v5_5_F2 在match_v5_4_F2的基础上
%               将最后“至少4个点”的步骤放到第4步删点之前，没有成功匹配的点不删，让之后的匹配可以匹配到之前没有用上的点。
% match_v7_1_F2 默认在历史版本目录的上一个版本的基础上修改
%               将第2步抛点和第3步检漏合并为双极线投票法，投票超过半数的点通过筛选。
% match_v7_2_F2 在投票过后，用3D重建点残差方法再抛点，以防几个镜头的射线处于同一极平面的情况。
% match_v7_3_F2 将第3步抛点放到第4步删点之后。
% match_v7_4_F2 在第1步选取第二个基准点pb时，除了保证同一个镜头中没有与pa匹配的其他点外，
%               增加另一个条件：pb是所有与pa匹配的点中匹配程度最高的，即点到极线的距离最小的。
%               计算pb的函数 onepointincam 改为 selectpb，同时实现匹配程度最高和同一镜头没有其他点的功能。
% match_v7_5_F2 在第1步选取第二个基准点pb时，添加处理同一个镜头中存在两个及以上匹配点的情况，
%               比较点到极线距离最小的和第二小的两个值的比值，
%               若比值小于0.6，则说明距离最小的点相对于第二小的点，匹配程度显著更高，选择距离最小的点作为匹配点；
%               若比值大于0.6，则说明距离最小的点相对于第二小的点，匹配程度没有显著的区别，无法判断此镜头中哪个点更好，都扔掉。
% match_v7_6_F2 在第1步选取第二个基准点pb时，pb的搜索条件由匹配程度最高，改为pb的镜头所在空间位置离pa镜头最近。
%               讲道理，由于物体遮挡的原因，空间位置靠近的镜头，同时看到一个marker点的概率大，
%               而空间位置相对距离较大（处于对角）的两个镜头，同时看到一个marker点的概率较小。
%               因此，在搜索与pa匹配的基准点时，考虑从pa所在镜头附近的镜头开始寻找，是合理的。
%               ！需要在selectpb中手动输入镜头间的相对位置信息，类似于链表的结构！
% match_v7_7_F2 在第1步选取第二个基准点pb时，处理同一个镜头存在两个及以上匹配点情况中增加一步：
%               1.用各个候选pb与pa构成基准点，与其他镜头进行匹配，若匹配出镜头最多的pb比第二多的多匹配到两个镜头及以上，
%                 则使用匹配出镜头最多的pb；若不满足上述条件，则进行第二步。
%               2.比较点到极线距离最小的和第二小的两个值的比值。
% match_v7_8_F2 在第1步选取第二个基准点pb时，若左右搜寻到了可匹配点（同一个镜头中没有其他点的），进行判断：
%               反过来将可匹配点与pa所在镜头匹配，若pa所在镜头中也只有pa和可匹配点匹配，则将这个可匹配点作为基准点；
%               若若pa所在镜头中除了pa外，还有其他的点与可匹配点匹配，则认为此可匹配点不够好，左右搜寻下一个镜头。
% match_v7_9_F2 用外参的信息计算镜头的相对空间位置，不再需要输入location_cam .
% match_v7_10_F2 修复了v7_9中左右搜索第二个匹配点步骤无效的问题.
%                投票时，将自己给自己投0票改为自己给自己投1票.
% match_v7_11_F2 点残差部分，在匹配组中点较少时，随机抽取基准点的意义不大。
%                但点少时，错点对整个匹配组重建效果的影响很大，需要针对性地处理点少（<=5）的情形。
%                点少时，整个匹配组中的点两两重建遍历。再聚类（？）计算出基准3D坐标和基准点。
% match_v8_0_F2  将删点步骤放到点残差抛点之后……
% match_v8_0_F3  函数名后缀带F3的，为点到极线的距离矩阵只计算上三角，提速……
% match_v9_0_F3  大改selectb，不再只寻找一个点作为selectb，能与pa共同匹配到很多点的点，都作为selectb进入基准点组.
% match_v9_1_F3  vote中添加判断一个镜头中匹配到两个及以上点的情况
% match_v9_2_F3  selectb中投票的阈值放宽
% match_v9_3_F3  改变了check_point中搜寻最小距离的方法
% match_v9_4_F3  2016.8.30
%                按算法实现文档中的顺序修改注释整理，
%                优化筛选函数screen的计算过程，
%                npoint中多添加一行“前面镜头点数和”信息，简化了cp2line函数中的运算，
%                增添参数camline，对应threshold_pass中每个点的镜头号和点号，第一行为镜头号，第二行为点号.
%                用camline上的操作替代line2cp函数的功能，
%                抛点check_point函数中，优化两两重建的方式，
%                优化寻找匹配组一个镜头中存在多个点的操作，删除和寻找分别在函数 DeleteMultiple 和 Multiple 中实现.
%                将第四步点残差抛点放到了匹配循环外，即先匹配出所有的匹配组，再抛点，将点残差抛点与匹配过程分离，可以用于分步式计算？
% match_v9_5_F3  2016.8.31
%                抛点check_points函数中，两两重建3D坐标后，不再用3D坐标间两两计算距离，
%                直接取xyz三个方向上的中位数作为标准点的3D坐标，计算重建的3D点到标准点的距离，
%                距离超过阈值的组合扔掉，剩下的3D坐标取平均值作为这一匹配组的3D坐标返回，
%                省去了之后计算3D坐标的时间.
%                抛点函数实际上的作用改为计算重建3D坐标……
% match_v9_6_F3  2016.9.2
%                试试看利用上一帧的2D坐标匹配呢……
% match_v10_0_F2 尝试使用交叉匹配方法减小计算量
% match_v10_1_F3 在得到与主点pa匹配较好的点构成的匹配组后，重建，将重建坐标投影到未匹配的镜头上寻找匹配点
% match_v10_3_F3 添加判断，重建的点分别来自不同镜头，即SelectPb的返回值中不包含匹配到多个点的镜头
% match_v10_5_F3 修改对同一镜头匹配到多个点的处理，将剔除改为比较
% match_v10_6_F3 重投影步骤(SelectOther)添加判断同一镜头匹配到多个点的情况
% match_v11_0_F3 2016.10.21
%                重建投影的方法在点数增多的情形下效果不如循环投票方法，核心算法改回投票法。
%                在 match_v9_4_F3 版本的算法基础上，沿用 match_v10_6_F3 中的数据结构。
% 


threshold2D = threshold ; % 点到极线距离的2D阈值
camN = length(uv) ;

npoint = zeros(1,camN);   % 每个镜头中拍到的点数，如npoint(3)的值为6，则表示3号镜头所拍到的点有6个。
sumnp = zeros(1,camN+1) ; % 第i个元素表示前i个镜头的点数和，不包含第i个镜头.
for i = 1:camN
    npoint(i) = length(uv{i})/2 ;
    sumnp(i+1) = sum(npoint(1:i)) ;
end

camline = zeros(2,sumnp(camN+1)) ;  % 对应threshold_pass中每个点的镜头号和点号，第一行为镜头号，第二行为点号
for i = 1:camN
    camline(:,sumnp(i)+1:sumnp(i+1)) = [i*ones(1,npoint(i)); 1:npoint(i)] ;
end

% 将uv变为2*n的形式
tuv = zeros(2,sumnp(camN+1)) ;
for i = 1:length(uv)
    tt = uv{i} ;
    tt = tt' ;
    tt = reshape(tt,2,length(tt)/2) ;
    tuv(:,sumnp(i)+1:sumnp(i+1)) = tt ;
end
uv = tuv ;

match_pass.i = 0 ;
match_pass.camN = camN ;
match_pass.npoint = npoint(:) ;
match_pass.sumnp = sumnp ;
match_pass.camline = camline ;
match_pass.uv = uv ;
pass = zeros(1,sumnp(camN+1)) ;

% ==================== Step.1 用点到极线的距离进行筛选 ====================
% threshold_pass N_allpoint*N_allpoint 所有镜头所有点之间的双极线筛选结果，通过阈值的为1，没有通过阈值的为0

threshold_pass0 = screen(uv,F,threshold2D,match_pass) ;
threshold_pass = threshold_pass0  ;


kpass = 0 ; %index_match_pass 匹配出来的组数
for ipa = 1:sumnp(camN+1) % index_point_a 第一个基准点在threshold_pass中的列标

    % =========================== Step.2 循环投票 ============================
    % 以pointa和pointb两个不同镜头中的相互匹配点为基准点，与其他镜头的点进行匹配；
    % 同时被pointa和pointb匹配到的点视为正确的匹配点。
    %（不同镜头中的两点，用双外极线法计算的z值通过阈值视为“被匹配到的两点”）
    % 同一个镜头中有两个点被匹配到的情况，怎么处理呢…… = =
    
    % ================ Step.2.0 选择基准点 =================
    pa = threshold_pass(ipa,:) ;
    if isnan(pa), continue; end
    
    pa1 = find(pa==1); %与pa匹配的点的下标
    ipb = selectpb(ipa,pa1,threshold_pass,match_pass) ;
    if isempty(ipb), continue; end
    mpass = [ipa, ipb] ;

    % ============ Step.2.1 ~2.2 组内全局循环投票 ============
    lengthpass = 0 ;
    k = 1 ; % 循环次数
    while lengthpass ~= length(mpass)
        lengthpass = length(mpass) ;
        mpass = vote(mpass,threshold_pass,match_pass) ;
        k = k + 1 ;
        if k > 3, break; end %防止投出去又投回来，反反复复的情况……正常情况下循环3~5次足够找到正确匹配组了
    end
    
    % 至少用mincam个镜头
    if length(mpass) < mincam
        continue;
    end
    
    % ============================ Step.3 删点 ==============================
    kpass = kpass + 1 ;
    pass(mpass) = kpass ;
    threshold_pass(mpass,:) = nan ;
    threshold_pass(:,mpass) = nan ;
    
end % ipa = 1:size(threshold_pass0,1)

match_pass.pass = pass ;


end %function match end

function threshold_pass = screen(uv,F,threshold2D,match_pass)
% 筛选
% 按threshold_pass的数据结构存储通过阈值点的镜头号、点号、z值（z=[uva 1]*F*[uvb 1]'）
% input:
% uv 1*camN cell型，单帧的各镜头uv坐标
% M n*11 各镜头的M矩阵
% threshold 1*1 点到极线距离的阈值
% output:
% Z N_allpoint*N_allpoint  所有镜头所有点之间的点到极线的距离
% threshold_pass N_allpoint*N_allpoint 所有镜头所有点之间的点到极线距离筛选结果，通过阈值的为1，没有通过阈值的为0

npoint = match_pass.npoint ;
sumnp = match_pass.sumnp ;
camN = match_pass.camN ;

Z = nan(sumnp(camN+1));

for icam = 1:camN-1
    iN = npoint(icam) ; %i号镜头拍到的点的个数
    if iN<1, continue; end
  
    for jcam = icam+1:camN
        jN = npoint(jcam) ; %j号镜头拍到的点的个数
        
        tF = F(icam*3-2:icam*3, jcam*3-2:jcam*3) ; %temp_F
        uvi = uv(:,sumnp(icam)+1:sumnp(icam+1)) ;
        uvj = uv(:,sumnp(jcam)+1:sumnp(jcam+1)) ;

        L = tF*[uvj; ones(1,jN)] ;
        tempF = [uvi', ones(iN,1)] * L ;
        L = sqrt(L(1,:).^2+L(2,:).^2) ;
        tempF = bsxfun( @rdivide, tempF, L ) ;  % tempF = tempF./repmat(L,iN,1) ;
        Z(sumnp(jcam)+1:sumnp(jcam+1), sumnp(icam)+1:sumnp(icam+1)) = tempF' ;
        Z(sumnp(icam)+1:sumnp(icam+1), sumnp(jcam)+1:sumnp(jcam+1)) = tempF ;
    end
end

Z = abs(Z) ;
threshold_pass = Z< threshold2D ;
threshold_pass = threshold_pass + eye(sumnp(camN+1)) ;
end % screen end

function ipb = selectpb(ipa,p,threshold_pass,match_pass)
% 查找p中是否存在能与pa一起匹配出较多点的点
% 若不存在，返回[].
% input:
% ipa    1*1       第一个基准点pa在大矩阵中的下标
% p      1行       与pa匹配的点在大矩阵中的下标
% threshold_pass    N_allpoint*N_allpoint 所有镜头所有点之间的点到极线距离筛选结果，通过阈值的为1，没有通过阈值的为0
% npoint 1行       每个镜头中拍到的点数，如npoint(3)的值为6，则表示3号镜头所拍到的点有6个。
% inexI  camN*10   镜头的内外参数
% output:
% ipb    1行       第二批基准点

sumnp = match_pass.sumnp ;
camline = match_pass.camline ;

if length(p)<2, ipb=[]; return; end  %如果只有ipa自己，就不用选了……

% 每一个备选点都和pa一起匹配，看哪一个匹配与pa一起匹配出来的点多
pp = zeros(1,length(p)); %备选点与pa一起匹配出的点数
for id = 1:length(p)
    if p(id) == ipa, continue; end
    tline =  nansum( threshold_pass([ipa,p(id)],:), 1 ) > 1  ; %得2票的点
    tcp = DeleteRepeated(camline(:,tline)) ;
    pp(id) = size(tcp,2) ;
end

temp = p(pp>=max(pp)*2/3) ; %选一些票数高的点……
tcp = DeleteRepeated(camline(:,temp)) ;
ipb = sumnp(tcp(1,:)) + tcp(2,:) ;
% ipb = cp2line(tcp,npoint) ;
end %selectpb

function cp = DeleteRepeated(cp)
% 删除一个镜头中存在多个点的镜头
% input
% cp        2行  一个匹配组的点按 [镜头号; 点号] 方式存储的形式
% output
% cp        2行  cam_point 2行 [镜头号; 点号] 

b = sort(cp(1,:)) ;
db = diff(b);
t = db~=0 ;
t = [true,t] & [t, true] ;
cp = cp(:,t) ;
% linep = cp2line(cp,npoint) ;

end

function b = Repeated(cp)
%提取存在两个及以上点的镜头号
% input
% cp        2行  一个匹配组的点按 [镜头号; 点号] 方式存储的形式
% output
% b         1行  cam# 存在多点的镜头号 

b = sort(cp(1,:));
db = diff(b);
b = b(db==0) ; %重复出现的元素
if isempty(b)
    return ;
end
db = diff(b);
d = db ~= 0;
d(numel(b)) = true; % Final element is always a member of unique list.
b = b(d);

end

function mpass = vote(mpass,threshold_pass,match_pass)
% 投票法筛选（点到极线的距离）

sumnp = match_pass.sumnp ;
camline = match_pass.camline ;

% 用投票的方法将组外得票多的点纳入匹配组中
% 已经在组中的点不被投票
tempthrepass = threshold_pass ; %temp_threshold_pass
for i = 1:length(mpass)
    icam = camline(1,mpass(i)) ;
    tempthrepass( :, sumnp(icam)+1: sumnp(icam+1) ) = nan ;
end

votenum = nansum( tempthrepass(mpass,:) ) ;

mpass = [mpass, find( votenum >= length(mpass)/2 ) ] ; %得票多于半数的点纳入组中
mpass = unique(mpass) ;


% 组内再投一次，将组内得票少的点剔除匹配组外
% 只用镜头中只有一个点的镜头投票，一个镜头多点的不参与
temp = DeleteRepeated(camline(:,mpass)) ;
temp = sumnp(temp(1,:)) + temp(2,:) ;

temp = threshold_pass(temp,mpass) ;
votenum = nansum(temp) ;

mpass = mpass( votenum > size(temp,1)/2 )  ; %得票多于半数的点留在组中


% 判断一个镜头中匹配到两个及以上点的情况
camlinepass = camline(:,mpass) ; 
iqmpi = Repeated(camlinepass) ; %存在两个及以上点的镜头号
for it = 1:length(iqmpi)
    voteit = votenum(:,camlinepass(1,:)==iqmpi(it)) ; %同一镜头中各点的得票数
    tempplace =  camlinepass(:,camlinepass(1,:)==iqmpi(it)) ; %同一镜头中各点的镜头号和点号
    [maxvote,maxplace] = max(voteit) ; 
    temp = voteit ;
    temp(temp==maxvote) = [] ;
    
    camlinepass(:,camlinepass(1,:)==iqmpi(it))=[] ; %将同一镜头中的各点从匹配组中删除，若存在显著更优的点再纳入回来
    if max(temp)/maxvote < 0.6 %得票第二多的点与最多的点票数做比较，若比例值小于0.6，则认为最多票数的点比其他点显著更优
        camlinepass(:,size(camlinepass,2)+1) = tempplace(:,maxplace) ;
    end
        
    mpass = sumnp(camlinepass(1,:)) + camlinepass(2,:) ;
end


end %vote




