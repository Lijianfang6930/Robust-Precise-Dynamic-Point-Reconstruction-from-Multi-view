function match_pass = match_v11_0_F3(uv,F,threshold,mincam)
% ����������ֵ����ƥ��
% ����Ϊ��
% 1.��ƥ�䣨˫���ߣ�
% 2.ͶƱ��ɸѡ��˫���ߣ�
% 3.���ɾ��
% 4.�׵㣨��в
%
% ���룺
% uv             1*camN    cell�ͣ���֡�����о�ͷuv��Ϣ��camNΪ��ͷ��������ͬ��
% M              camN*11   ÿ����ͷ��������M����
% threshold      1*1       �㵽���߾������ֵ
% mincam         ������ֵ   ƥ��������ʹ�õľ�ͷ����
% inexI          camN*10   ��ͷ���������
% �����
% match_pass 1*n      cell�ͣ�ÿ��Ԫ�ر�ʾƥ�����һ����
%                     cell�еĸ�ʽ��2*N����һ��Ϊ��ͷ�ţ��ڶ���Ϊ��ͷ�еĵ��

% ��ʷ�汾��
% match_v1_1    ����������ȥ���ڴ��Ĳ��裬�Ƚ������ؽ��Ĳ��죬ȥ���ؽ�Ч�����õ�uv�㣬��������ǵ�в����
% match_v1_2    �׵�ʹ���߲в�
% match_v1_3    ����match_v1_1�ĵ�в����в�check_points��step2�У�
%               �����ؽ�Ѱ�һ�׼uv��ʱ��Ѱ��˳���Ϊ���򣨸о���û��ʲô�á�����.
% match_v1_4    ��д��check_points���޸��˵�в�ıȽϷ�ʽ.
% match_v1_5    ��match_v1_4�Ļ����ϣ��޸�����ƥ�������Ϊ��
%               1.��ƥ�䣨˫���ߣ�
%               2.�����׵㣨��в
%               3.��©����в
%               4.���ɾ��
%               �ؽ��������������ˣ����ǵ㶶������
%               �����ԭ������ֵѡ���ˣ���С��ֵ�����þ�û��ô������
% match_v5      ��match_v1_5�Ļ�����
%               ��threshold_pass�����ݽṹ��ɴ�������ʽ.
%               �к��е�Ԫ�ؾ��ǣ� C1P1 C1P2 ... C1Pn_1 C2P1 C2P2 ... C2Pn_2 C3P1 ... CnP1 CnP2 ... CnPn_n
%               ����CiPj��ʾ��i�ž�ͷ�ĵ�j���㣬n_k��ʾ��k����ͷ�е������Ϊn_k.
%               threshold_pass�ĵ�C1P1�С�C3P1�б�ʾ��1�ž�ͷ��1�ŵ����3�ž�ͷ��1�ŵ�֮���˫���߷�������ֵz.
%               �Ӻ�����Ӧ�����˸ı䣬��ͬʱƥ�䵽ͬһ����ͷ������������û��������
% match_v5_1_F2 ��match_v5�Ļ�����
%               ɸѡ����screen�У���Ϊ�жϵ㵽���ߵľ�����ֵ��ͳһ��ֵΪ�̶�ֵ������ÿ������ͷ���㲻ͬ����ֵ
%               ��������׺��F2�ģ���Ϊ�㵽���ߵľ�����ֵ
% match_v5_2_F2 ��match_v5_1_F2�Ļ�����
%               ˫���߷�����ƥ��Ĳ����У�Ӧ��֤��׼��pb���ڵľ�ͷ�У�û��������paƥ��ĵ�
% match_v5_4_F2 ��match_v5_2_F2�Ļ�����
%               match_v5_2_F2�е�3����©�Ŀ�ʼ���֣���ɾ����ƥ��ĵ㣻
%               ��match_v5_4_F2�У���ɾ����Ĳ���ŵ���2���Ŀ�ʼ���֣���ֹ��3���������һ�.
% match_v5_5_F2 ��match_v5_4_F2�Ļ�����
%               ���������4���㡱�Ĳ���ŵ���4��ɾ��֮ǰ��û�гɹ�ƥ��ĵ㲻ɾ����֮���ƥ�����ƥ�䵽֮ǰû�����ϵĵ㡣
% match_v7_1_F2 Ĭ������ʷ�汾Ŀ¼����һ���汾�Ļ������޸�
%               ����2���׵�͵�3����©�ϲ�Ϊ˫����ͶƱ����ͶƱ���������ĵ�ͨ��ɸѡ��
% match_v7_2_F2 ��ͶƱ������3D�ؽ���в�����׵㣬�Է�������ͷ�����ߴ���ͬһ��ƽ��������
% match_v7_3_F2 ����3���׵�ŵ���4��ɾ��֮��
% match_v7_4_F2 �ڵ�1��ѡȡ�ڶ�����׼��pbʱ�����˱�֤ͬһ����ͷ��û����paƥ����������⣬
%               ������һ��������pb��������paƥ��ĵ���ƥ��̶���ߵģ����㵽���ߵľ�����С�ġ�
%               ����pb�ĺ��� onepointincam ��Ϊ selectpb��ͬʱʵ��ƥ��̶���ߺ�ͬһ��ͷû��������Ĺ��ܡ�
% match_v7_5_F2 �ڵ�1��ѡȡ�ڶ�����׼��pbʱ����Ӵ���ͬһ����ͷ�д�������������ƥ���������
%               �Ƚϵ㵽���߾�����С�ĺ͵ڶ�С������ֵ�ı�ֵ��
%               ����ֵС��0.6����˵��������С�ĵ�����ڵڶ�С�ĵ㣬ƥ��̶��������ߣ�ѡ�������С�ĵ���Ϊƥ��㣻
%               ����ֵ����0.6����˵��������С�ĵ�����ڵڶ�С�ĵ㣬ƥ��̶�û�������������޷��жϴ˾�ͷ���ĸ�����ã����ӵ���
% match_v7_6_F2 �ڵ�1��ѡȡ�ڶ�����׼��pbʱ��pb������������ƥ��̶���ߣ���Ϊpb�ľ�ͷ���ڿռ�λ����pa��ͷ�����
%               ���������������ڵ���ԭ�򣬿ռ�λ�ÿ����ľ�ͷ��ͬʱ����һ��marker��ĸ��ʴ�
%               ���ռ�λ����Ծ���ϴ󣨴��ڶԽǣ���������ͷ��ͬʱ����һ��marker��ĸ��ʽ�С��
%               ��ˣ���������paƥ��Ļ�׼��ʱ�����Ǵ�pa���ھ�ͷ�����ľ�ͷ��ʼѰ�ң��Ǻ���ġ�
%               ����Ҫ��selectpb���ֶ����뾵ͷ������λ����Ϣ������������Ľṹ��
% match_v7_7_F2 �ڵ�1��ѡȡ�ڶ�����׼��pbʱ������ͬһ����ͷ��������������ƥ������������һ����
%               1.�ø�����ѡpb��pa���ɻ�׼�㣬��������ͷ����ƥ�䣬��ƥ�����ͷ����pb�ȵڶ���Ķ�ƥ�䵽������ͷ�����ϣ�
%                 ��ʹ��ƥ�����ͷ����pb������������������������еڶ�����
%               2.�Ƚϵ㵽���߾�����С�ĺ͵ڶ�С������ֵ�ı�ֵ��
% match_v7_8_F2 �ڵ�1��ѡȡ�ڶ�����׼��pbʱ����������Ѱ���˿�ƥ��㣨ͬһ����ͷ��û��������ģ��������жϣ�
%               ����������ƥ�����pa���ھ�ͷƥ�䣬��pa���ھ�ͷ��Ҳֻ��pa�Ϳ�ƥ���ƥ�䣬�������ƥ�����Ϊ��׼�㣻
%               ����pa���ھ�ͷ�г���pa�⣬���������ĵ����ƥ���ƥ�䣬����Ϊ�˿�ƥ��㲻���ã�������Ѱ��һ����ͷ��
% match_v7_9_F2 ����ε���Ϣ���㾵ͷ����Կռ�λ�ã�������Ҫ����location_cam .
% match_v7_10_F2 �޸���v7_9�����������ڶ���ƥ��㲽����Ч������.
%                ͶƱʱ�����Լ����Լ�Ͷ0Ʊ��Ϊ�Լ����Լ�Ͷ1Ʊ.
% match_v7_11_F2 ��в�֣���ƥ�����е����ʱ�������ȡ��׼������岻��
%                ������ʱ����������ƥ�����ؽ�Ч����Ӱ��ܴ���Ҫ����Եش�����٣�<=5�������Ρ�
%                ����ʱ������ƥ�����еĵ������ؽ��������پ��ࣨ�����������׼3D����ͻ�׼�㡣
% match_v8_0_F2  ��ɾ�㲽��ŵ���в��׵�֮�󡭡�
% match_v8_0_F3  ��������׺��F3�ģ�Ϊ�㵽���ߵľ������ֻ���������ǣ����١���
% match_v9_0_F3  ���selectb������ֻѰ��һ������Ϊselectb������pa��ͬƥ�䵽�ܶ��ĵ㣬����Ϊselectb�����׼����.
% match_v9_1_F3  vote������ж�һ����ͷ��ƥ�䵽���������ϵ�����
% match_v9_2_F3  selectb��ͶƱ����ֵ�ſ�
% match_v9_3_F3  �ı���check_point����Ѱ��С����ķ���
% match_v9_4_F3  2016.8.30
%                ���㷨ʵ���ĵ��е�˳���޸�ע������
%                �Ż�ɸѡ����screen�ļ�����̣�
%                npoint�ж����һ�С�ǰ�澵ͷ�����͡���Ϣ������cp2line�����е����㣬
%                �������camline����Ӧthreshold_pass��ÿ����ľ�ͷ�ź͵�ţ���һ��Ϊ��ͷ�ţ��ڶ���Ϊ���.
%                ��camline�ϵĲ������line2cp�����Ĺ��ܣ�
%                �׵�check_point�����У��Ż������ؽ��ķ�ʽ��
%                �Ż�Ѱ��ƥ����һ����ͷ�д��ڶ����Ĳ�����ɾ����Ѱ�ҷֱ��ں��� DeleteMultiple �� Multiple ��ʵ��.
%                �����Ĳ���в��׵�ŵ���ƥ��ѭ���⣬����ƥ������е�ƥ���飬���׵㣬����в��׵���ƥ����̷��룬�������ڷֲ�ʽ���㣿
% match_v9_5_F3  2016.8.31
%                �׵�check_points�����У������ؽ�3D����󣬲�����3D���������������룬
%                ֱ��ȡxyz���������ϵ���λ����Ϊ��׼���3D���꣬�����ؽ���3D�㵽��׼��ľ��룬
%                ���볬����ֵ������ӵ���ʣ�µ�3D����ȡƽ��ֵ��Ϊ��һƥ�����3D���귵�أ�
%                ʡȥ��֮�����3D�����ʱ��.
%                �׵㺯��ʵ���ϵ����ø�Ϊ�����ؽ�3D���ꡭ��
% match_v9_6_F3  2016.9.2
%                ���Կ�������һ֡��2D����ƥ���ء���
% match_v10_0_F2 ����ʹ�ý���ƥ�䷽����С������
% match_v10_1_F3 �ڵõ�������paƥ��Ϻõĵ㹹�ɵ�ƥ������ؽ������ؽ�����ͶӰ��δƥ��ľ�ͷ��Ѱ��ƥ���
% match_v10_3_F3 ����жϣ��ؽ��ĵ�ֱ����Բ�ͬ��ͷ����SelectPb�ķ���ֵ�в�����ƥ�䵽�����ľ�ͷ
% match_v10_5_F3 �޸Ķ�ͬһ��ͷƥ�䵽�����Ĵ������޳���Ϊ�Ƚ�
% match_v10_6_F3 ��ͶӰ����(SelectOther)����ж�ͬһ��ͷƥ�䵽���������
% match_v11_0_F3 2016.10.21
%                �ؽ�ͶӰ�ķ����ڵ��������������Ч������ѭ��ͶƱ�����������㷨�Ļ�ͶƱ����
%                �� match_v9_4_F3 �汾���㷨�����ϣ����� match_v10_6_F3 �е����ݽṹ��
% 


threshold2D = threshold ; % �㵽���߾����2D��ֵ
camN = length(uv) ;

npoint = zeros(1,camN);   % ÿ����ͷ���ĵ��ĵ�������npoint(3)��ֵΪ6�����ʾ3�ž�ͷ���ĵ��ĵ���6����
sumnp = zeros(1,camN+1) ; % ��i��Ԫ�ر�ʾǰi����ͷ�ĵ����ͣ���������i����ͷ.
for i = 1:camN
    npoint(i) = length(uv{i})/2 ;
    sumnp(i+1) = sum(npoint(1:i)) ;
end

camline = zeros(2,sumnp(camN+1)) ;  % ��Ӧthreshold_pass��ÿ����ľ�ͷ�ź͵�ţ���һ��Ϊ��ͷ�ţ��ڶ���Ϊ���
for i = 1:camN
    camline(:,sumnp(i)+1:sumnp(i+1)) = [i*ones(1,npoint(i)); 1:npoint(i)] ;
end

% ��uv��Ϊ2*n����ʽ
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

% ==================== Step.1 �õ㵽���ߵľ������ɸѡ ====================
% threshold_pass N_allpoint*N_allpoint ���о�ͷ���е�֮���˫����ɸѡ�����ͨ����ֵ��Ϊ1��û��ͨ����ֵ��Ϊ0

threshold_pass0 = screen(uv,F,threshold2D,match_pass) ;
threshold_pass = threshold_pass0  ;


kpass = 0 ; %index_match_pass ƥ�����������
for ipa = 1:sumnp(camN+1) % index_point_a ��һ����׼����threshold_pass�е��б�

    % =========================== Step.2 ѭ��ͶƱ ============================
    % ��pointa��pointb������ͬ��ͷ�е��໥ƥ���Ϊ��׼�㣬��������ͷ�ĵ����ƥ�䣻
    % ͬʱ��pointa��pointbƥ�䵽�ĵ���Ϊ��ȷ��ƥ��㡣
    %����ͬ��ͷ�е����㣬��˫�⼫�߷������zֵͨ����ֵ��Ϊ����ƥ�䵽�����㡱��
    % ͬһ����ͷ���������㱻ƥ�䵽���������ô�����ء��� = =
    
    % ================ Step.2.0 ѡ���׼�� =================
    pa = threshold_pass(ipa,:) ;
    if isnan(pa), continue; end
    
    pa1 = find(pa==1); %��paƥ��ĵ���±�
    ipb = selectpb(ipa,pa1,threshold_pass,match_pass) ;
    if isempty(ipb), continue; end
    mpass = [ipa, ipb] ;

    % ============ Step.2.1 ~2.2 ����ȫ��ѭ��ͶƱ ============
    lengthpass = 0 ;
    k = 1 ; % ѭ������
    while lengthpass ~= length(mpass)
        lengthpass = length(mpass) ;
        mpass = vote(mpass,threshold_pass,match_pass) ;
        k = k + 1 ;
        if k > 3, break; end %��ֹͶ��ȥ��Ͷ��������������������������������ѭ��3~5���㹻�ҵ���ȷƥ������
    end
    
    % ������mincam����ͷ
    if length(mpass) < mincam
        continue;
    end
    
    % ============================ Step.3 ɾ�� ==============================
    kpass = kpass + 1 ;
    pass(mpass) = kpass ;
    threshold_pass(mpass,:) = nan ;
    threshold_pass(:,mpass) = nan ;
    
end % ipa = 1:size(threshold_pass0,1)

match_pass.pass = pass ;


end %function match end

function threshold_pass = screen(uv,F,threshold2D,match_pass)
% ɸѡ
% ��threshold_pass�����ݽṹ�洢ͨ����ֵ��ľ�ͷ�š���š�zֵ��z=[uva 1]*F*[uvb 1]'��
% input:
% uv 1*camN cell�ͣ���֡�ĸ���ͷuv����
% M n*11 ����ͷ��M����
% threshold 1*1 �㵽���߾������ֵ
% output:
% Z N_allpoint*N_allpoint  ���о�ͷ���е�֮��ĵ㵽���ߵľ���
% threshold_pass N_allpoint*N_allpoint ���о�ͷ���е�֮��ĵ㵽���߾���ɸѡ�����ͨ����ֵ��Ϊ1��û��ͨ����ֵ��Ϊ0

npoint = match_pass.npoint ;
sumnp = match_pass.sumnp ;
camN = match_pass.camN ;

Z = nan(sumnp(camN+1));

for icam = 1:camN-1
    iN = npoint(icam) ; %i�ž�ͷ�ĵ��ĵ�ĸ���
    if iN<1, continue; end
  
    for jcam = icam+1:camN
        jN = npoint(jcam) ; %j�ž�ͷ�ĵ��ĵ�ĸ���
        
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
% ����p���Ƿ��������paһ��ƥ����϶��ĵ�
% �������ڣ�����[].
% input:
% ipa    1*1       ��һ����׼��pa�ڴ�����е��±�
% p      1��       ��paƥ��ĵ��ڴ�����е��±�
% threshold_pass    N_allpoint*N_allpoint ���о�ͷ���е�֮��ĵ㵽���߾���ɸѡ�����ͨ����ֵ��Ϊ1��û��ͨ����ֵ��Ϊ0
% npoint 1��       ÿ����ͷ���ĵ��ĵ�������npoint(3)��ֵΪ6�����ʾ3�ž�ͷ���ĵ��ĵ���6����
% inexI  camN*10   ��ͷ���������
% output:
% ipb    1��       �ڶ�����׼��

sumnp = match_pass.sumnp ;
camline = match_pass.camline ;

if length(p)<2, ipb=[]; return; end  %���ֻ��ipa�Լ����Ͳ���ѡ�ˡ���

% ÿһ����ѡ�㶼��paһ��ƥ�䣬����һ��ƥ����paһ��ƥ������ĵ��
pp = zeros(1,length(p)); %��ѡ����paһ��ƥ����ĵ���
for id = 1:length(p)
    if p(id) == ipa, continue; end
    tline =  nansum( threshold_pass([ipa,p(id)],:), 1 ) > 1  ; %��2Ʊ�ĵ�
    tcp = DeleteRepeated(camline(:,tline)) ;
    pp(id) = size(tcp,2) ;
end

temp = p(pp>=max(pp)*2/3) ; %ѡһЩƱ���ߵĵ㡭��
tcp = DeleteRepeated(camline(:,temp)) ;
ipb = sumnp(tcp(1,:)) + tcp(2,:) ;
% ipb = cp2line(tcp,npoint) ;
end %selectpb

function cp = DeleteRepeated(cp)
% ɾ��һ����ͷ�д��ڶ����ľ�ͷ
% input
% cp        2��  һ��ƥ����ĵ㰴 [��ͷ��; ���] ��ʽ�洢����ʽ
% output
% cp        2��  cam_point 2�� [��ͷ��; ���] 

b = sort(cp(1,:)) ;
db = diff(b);
t = db~=0 ;
t = [true,t] & [t, true] ;
cp = cp(:,t) ;
% linep = cp2line(cp,npoint) ;

end

function b = Repeated(cp)
%��ȡ�������������ϵ�ľ�ͷ��
% input
% cp        2��  һ��ƥ����ĵ㰴 [��ͷ��; ���] ��ʽ�洢����ʽ
% output
% b         1��  cam# ���ڶ��ľ�ͷ�� 

b = sort(cp(1,:));
db = diff(b);
b = b(db==0) ; %�ظ����ֵ�Ԫ��
if isempty(b)
    return ;
end
db = diff(b);
d = db ~= 0;
d(numel(b)) = true; % Final element is always a member of unique list.
b = b(d);

end

function mpass = vote(mpass,threshold_pass,match_pass)
% ͶƱ��ɸѡ���㵽���ߵľ��룩

sumnp = match_pass.sumnp ;
camline = match_pass.camline ;

% ��ͶƱ�ķ����������Ʊ��ĵ�����ƥ������
% �Ѿ������еĵ㲻��ͶƱ
tempthrepass = threshold_pass ; %temp_threshold_pass
for i = 1:length(mpass)
    icam = camline(1,mpass(i)) ;
    tempthrepass( :, sumnp(icam)+1: sumnp(icam+1) ) = nan ;
end

votenum = nansum( tempthrepass(mpass,:) ) ;

mpass = [mpass, find( votenum >= length(mpass)/2 ) ] ; %��Ʊ���ڰ����ĵ���������
mpass = unique(mpass) ;


% ������Ͷһ�Σ������ڵ�Ʊ�ٵĵ��޳�ƥ������
% ֻ�þ�ͷ��ֻ��һ����ľ�ͷͶƱ��һ����ͷ���Ĳ�����
temp = DeleteRepeated(camline(:,mpass)) ;
temp = sumnp(temp(1,:)) + temp(2,:) ;

temp = threshold_pass(temp,mpass) ;
votenum = nansum(temp) ;

mpass = mpass( votenum > size(temp,1)/2 )  ; %��Ʊ���ڰ����ĵ���������


% �ж�һ����ͷ��ƥ�䵽���������ϵ�����
camlinepass = camline(:,mpass) ; 
iqmpi = Repeated(camlinepass) ; %�������������ϵ�ľ�ͷ��
for it = 1:length(iqmpi)
    voteit = votenum(:,camlinepass(1,:)==iqmpi(it)) ; %ͬһ��ͷ�и���ĵ�Ʊ��
    tempplace =  camlinepass(:,camlinepass(1,:)==iqmpi(it)) ; %ͬһ��ͷ�и���ľ�ͷ�ź͵��
    [maxvote,maxplace] = max(voteit) ; 
    temp = voteit ;
    temp(temp==maxvote) = [] ;
    
    camlinepass(:,camlinepass(1,:)==iqmpi(it))=[] ; %��ͬһ��ͷ�еĸ����ƥ������ɾ�����������������ŵĵ����������
    if max(temp)/maxvote < 0.6 %��Ʊ�ڶ���ĵ������ĵ�Ʊ�����Ƚϣ�������ֵС��0.6������Ϊ���Ʊ���ĵ����������������
        camlinepass(:,size(camlinepass,2)+1) = tempplace(:,maxplace) ;
    end
        
    mpass = sumnp(camlinepass(1,:)) + camlinepass(2,:) ;
end


end %vote




