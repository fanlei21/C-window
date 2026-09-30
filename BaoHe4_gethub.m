%% ========== 土地利用掩膜 + 分区统计（只保留指定地类） ==========
clear; clc; close all;

% ===================== 1. 【你只需要改这里】设置参数 =====================
% 你要保留的土地利用类型（改成你自己的！）
keep_land_classes = [40,50,60,61,62,70,71,72,80,81,82,90,100,110,120,121,122,130,140,150,151,152,153,160,170,180];  % 例如：1=耕地,2=森林,5=草地，按需修改

% 文件路径
land_path = 'H:\数据库\ESA_CCI_LandCoverData\ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7.tif';   % 土地利用掩膜
A_path    = 'H:\数据库\气候类型数据\koppen_geiger_tif\1991_2020\koppen_geiger_0p5.tif';          % 分类图 A
% B_path    = 'H:\数据库\CMIP6_SSP585_Soil_Carbon_data\TIFF_输出2\cSoil_mean_2094_change.tif';          % 你的数据图 B

% ======================================================================
% 1. 读取 A 图（低分辨率，作为目标网格）
[A, ~] = geotiffread(A_path);
[~, R_A] = geotiffread('H:\数据库\气候类型数据\koppen_geiger_tif\1991_2020\temp0.5.tif');
targetSize = size(A);  % 目标尺寸 = A图的大小
[rowsA, colsA] = size(A);
% 生成 A 图经纬度（后续给 B 图匹配用）
lat_A = linspace(R_A.LatitudeLimits(1), R_A.LatitudeLimits(2), rowsA);
lon_A = linspace(R_A.LongitudeLimits(1), R_A.LongitudeLimits(2), colsA);
[LON_A, LAT_A] = meshgrid(lon_A, lat_A);

% ========== 步骤2：读取 300m 土地图，直接粗化==========
fprintf('正在粗化 300m 土地图...\n');
[land, ~] = geotiffread(land_path);

% 直接缩放
land = imresize(land, targetSize, 'nearest');  % 最近邻插值

% 2. 生成掩膜：只保留指定土地类型
land_mask = ismember(land, keep_land_classes);

% ========== 步骤3：B 图按经纬度匹配到 A 图（正确空间对齐） ==========
% type={'','Fast','Medium','Slow'};
% yy=[2026,2039,2060,2078,2094];
type={'1.5C','2.0C','3.0C','4.0C','5.0C'};
yy={'','TOC'};
for tt=1:5 % tt=1:4
    for ii = 1:2 %ii = 1:5
        %         B_path = fullfile('H:\数据库\CMIP6_SSP585_Soil_Carbon_data\TIFF_输出2', ['cSoil',type{tt},sprintf('_mean_%d_change.tif', yy(ii))]);
        B_path = fullfile('H:\数据库\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4', ...
            ['全球碳库SSP预测30_',type{tt},'增量',yy{ii},'.tif']);
        fprintf('B 图空间匹配中...\n');
        %         ==================================================================
        %         [B, R_B] = geotiffread(B_path);
        %         [rowsB, colsB] = size(B);
        %
        %         % 1. B图原始经纬度（0~360）
        %         lon_B_0_360 = linspace(R_B.LongitudeLimits(1), R_B.LongitudeLimits(2), colsB);
        %         lat_B = linspace(R_B.LatitudeLimits(1), R_B.LatitudeLimits(2), rowsB);
        %
        %         % 2. 经度转换：0~360 → -180~180
        %         lon_B = lon_B_0_360;
        %         lon_B(lon_B > 180) = lon_B(lon_B > 180) - 360;
        %
        %         % 3. 关键：B矩阵+经度左右翻转，匹配-180~180顺序
        %         [lon_B, idx_sort] = sort(lon_B);
        %         B = B(:, idx_sort);
        %         [LON_B, LAT_B] = meshgrid(lon_B, lat_B);
        %===========================================================================
        
        % ========== 步骤3：读取 B 图（自适应2D/3D，-180~180标准经纬度）==========
        [B, R_B] = geotiffread(B_path);
        [rowsB, colsB, ~] = size(B);  % ~ 自动忽略第三维（时间）
        
        % 直接生成标准经纬度（ -180~180）
        lat_B = linspace(R_B.LatitudeLimits(1), R_B.LatitudeLimits(2), rowsB);
        lon_B = linspace(R_B.LongitudeLimits(1), R_B.LongitudeLimits(2), colsB);
        [LON_B, LAT_B] = meshgrid(lon_B, lat_B);
        
        % 插值到 A 图（自动支持 2D 和 3D 数据）
        if ndims(B) == 2
            % 2D 数据
            B_on_A = interp2(LON_B, LAT_B, double(B), LON_A, LAT_A, 'linear', NaN);
        else
            % 3D 数据（自动逐时间层插值）
            B_on_A = zeros(size(LON_A,1), size(LON_A,2), size(B,3));
            for t = 1:size(B,3)
                B_on_A(:,:,t) = interp2(LON_B, LAT_B, double(B(:,:,t)), LON_A, LAT_A, 'linear', NaN);
            end
        end
        
        % ========== 步骤4：掩膜 ==========
        A_masked = A;
        A_masked(~land_mask) = NaN;

        % 如果B是3D，掩膜也要逐层
        if ndims(B_on_A) == 2
            B_masked = B_on_A;
            B_masked(~land_mask) = NaN;
        else
            B_masked = B_on_A;
            for t = 1:size(B_masked,3)
                tmp = B_masked(:,:,t);
                tmp(~land_mask) = NaN;
                B_masked(:,:,t) = tmp;
            end
        end
        

        classes = unique(A_masked);
        classes(isnan(classes)) = [];
        n_classes = length(classes);
        n_layers = size(B_masked, 3); % 层数

        % 预先创建矩阵：[分类数 × 6指标 × 层数]
        result_table = zeros(n_classes, 6, n_layers);
        
        for t = 1:n_layers
            B_stat = B_masked(:,:,t);  % 取当前层
    
            for i = 1:n_classes
                cls = classes(i);
                vals = B_stat(A_masked == cls);
                vals = vals(~isnan(vals));
        
                if isempty(vals)
                    result_table(i,:,t) = NaN;
                    continue;
                end
        
                % 一行6个指标              
                result_table(i,1,t) = double(cls);             % 分类
                result_table(i,2,t) = double(length(vals));    % 像元数
                result_table(i,3,t) = mean(vals,'omitnan');    % 均值
                result_table(i,4,t) = std(vals,'omitnan');     % 标准差
                result_table(i,5,t) = prctile(vals,5);         % P05
                result_table(i,6,t) = prctile(vals,95);        % P95
            end
        end
        
        %         mat_name = fullfile('H:\数据库\分区统计', ['cSoil',type{tt},sprintf('_mean_%d_change.mat', yy(ii))]);
        mat_name = fullfile('H:\数据库\分区统计', ['全球碳库SSP预测30_',type{tt},'增量',yy{ii},'.mat']);
        save(mat_name, 'result_table', '-v6');
    end
end


% % ========== 输出结果 ==========
% fprintf('-----------------------------------------------------\n');
% fprintf('分类\t像元数\t均值\t标准差\tp05\tp95\n');
% fprintf('-----------------------------------------------------\n');
% disp(result_table);
% 
% % 保存到 Excel（带表头）
% result_save = array2table(result_table,...
%     'VariableNames',{'Class','Count','Mean','Std','P05','P95'});
% mat_name = fullfile('H:\数据库\分区统计', ['cSoil',type,sprintf('_mean_%d_change.xlsx', yy(ii))]);
% writetable(result_save, excel_name);
% 
% fprintf('? 已保存到 Excel！\n');

% ========== 绘制柱状图：均值±标准差 + 气候类型名称标签 ==========
figure('Color','w','Position',[100,100,1200,600]);

% 1. 气候分类映射
climate_code = {
    1,'Af','Tropical, rainforest'
    2,'Am','Tropical, monsoon'
    3,'Aw','Tropical, savannah'
    4,'BWh','Arid, desert, hot'
    5,'BWk','Arid, desert, cold'
    6,'BSh','Arid, steppe, hot'
    7,'BSk','Arid, steppe, cold'
    8,'Csa','Temperate, dry summer, hot summer'
    9,'Csb','Temperate, dry summer, warm summer'
    10,'Csc','Temperate, dry summer, cold summer'
    11,'Cwa','Temperate, dry winter, hot summer'
    12,'Cwb','Temperate, dry winter, warm summer'
    13,'Cwc','Temperate, dry winter, cold summer'
    14,'Cfa','Temperate, no dry season, hot summer'
    15,'Cfb','Temperate, no dry season, warm summer'
    16,'Cfc','Temperate, no dry season, cold summer'
    17,'Dsa','Cold, dry summer, hot summer'
    18,'Dsb','Cold, dry summer, warm summer'
    19,'Dsc','Cold, dry summer, cold summer'
    20,'Dsd','Cold, dry summer, very cold winter'
    21,'Dwa','Cold, dry winter, hot summer'
    22,'Dwb','Cold, dry winter, warm summer'
    23,'Dwc','Cold, dry winter, cold summer'
    24,'Dwd','Cold, dry winter, very cold winter'
    25,'Dfa','Cold, no dry season, hot summer'
    26,'Dfb','Cold, no dry season, warm summer'
    27,'Dfc','Cold, no dry season, cold summer'
    28,'Dfd','Cold, no dry season, very cold winter'
    29,'ET','Polar, tundra'
    30,'EF','Polar, frost'
};

% 2. 提取当前统计结果
has_nan = any(isnan(result_table), 2);   
result_table(has_nan, :) = [];           

cls_num = result_table(:,1);   % 分类编号
means   = result_table(:,3);   % 均值
stds    = result_table(:,4);   % 标准差
n_vals  = result_table(:,2);   % ? 这里是你要的 N（像元数，第2列）
n_cls   = length(cls_num);

% 3. 匹配简写标签
x_tick_label = cell(n_cls,1);
for i = 1:n_cls
    idx = cell2mat(climate_code(:,1)) == cls_num(i);
    x_tick_label{i} = climate_code{idx,2};
end

% 4. 绘制柱状图+误差棒
bar(1:n_cls, means, 0.7, 'FaceColor',[0.25 0.6 0.85],'EdgeColor','k');
hold on;
errorbar(1:n_cls, means, stds, 'k','LineStyle','none','LineWidth',1.5);

% ========== 关键：在每个柱子上方标注 N 值 ==========
for i = 1:n_cls
    % 标注位置：柱子正上方
    text(i, means(i)+max(stds)*0.1, ...
         sprintf('N=%d', n_vals(i)), ...
         'HorizontalAlignment','center',...
         'FontSize',10,'FontWeight','bold','Color','k');
end

% 5. 坐标轴与标签设置
xticks(1:n_cls);
xticklabels(x_tick_label);
xtickangle(45);   

xlabel('Kppen Climate Class','FontSize',13,'FontWeight','bold');
ylabel('Mean Value','FontSize',13,'FontWeight','bold');
title('Mean ± Standard Deviation of Soil Carbon by Climate Type','FontSize',15,'FontWeight','bold');
grid on; grid minor;
set(gca,'FontSize',11,'Box','on');


% 6. 保存高清图
saveas(gcf,'气候分区_均值标准差柱状图.tif');
saveas(gcf,'气候分区_均值标准差柱状图.png');
disp('?柱状图已保存！');


%% CMIP6 cSoil 碳库数据批量读取+模式统计+年平均+多模式集合平均
clear; clc; close all;

% 1. 设置文件夹路径
folder_path = 'H:\数据库\CMIP6_SSP585_Soil_Carbon_data';

% % 2. 获取所有 nc 文件
% fileList = dir(fullfile(folder_path,'cSoil_*.nc'));  % 这里加了筛选！
% N_file = length(fileList);
% fprintf('共找到 %d 个 cSoil NC 文件\n', N_file);

% 2. 获取所有 nc 文件 cSoilFast
fileList = dir(fullfile(folder_path,'cSoilMedium_*.nc'));  % 这里加了筛选！
N_file = length(fileList);
fprintf('共找到 %d 个 cSoilSlow NC 文件\n', N_file);

if N_file == 0
    error('未找到 cSoil 开头的 NC 文件，请检查路径或文件名');
end

% 3. 初始化存储
% 模式信息
model_list = {};       % 模式名称
model_info_table = []; % 模式属性汇总表
% 年尺度数据
year_all = [];         % 所有年份
cSoil_year_all = {};   % 每个模式的逐年数据 [year, cSoil_mean]
attr_all = [];         % 存储每个文件的属性

% 4. 循环读取每个文件
for i = 1:N_file
    fname = fileList(i).name;
    ncpath = fullfile(folder_path,fname);
    fprintf('正在读取：%s\n',fname);
    
    % ========== 4.1 解析文件名信息 ==========
    split_str = split(fname,'_');
    varname = split_str{1};       % cSoil
    table_id = split_str{2};      % Emon / Eyr
    %model = split_str{3};         % 模式名
    model = strrep(split_str{3}, '-', '_');  % 关键修复：- → _
    experiment = split_str{4};    % ssp585
    run = split_str{5};           % r1i1p1f1
    grid_type = split_str{6};          % gn/gr
    
    % 提取时间范围
    time_str = erase(split_str{7},'.nc');
    if strcmp(table_id,'Emon') || strcmp(table_id,'Lmon')
        % Emon格式：201501-210012
        t_start = str2double(time_str(1:4));
        t_end   = str2double(time_str(8:11));
        m_start = str2double(time_str(5:6));
        n_end   = str2double(time_str(12:13));
    elseif strcmp(table_id,'Eyr')
        % Eyr格式：2015-2100
        t_start = str2double(extractBefore(time_str,'-'));
        t_end   = str2double(extractAfter(time_str,'-'));
    end
    
    % 存入模式列表
    if ~ismember(model,model_list)
        model_list{end+1} = model;
    end
    
    % ========== 4.2 读取 cSoil 变量属性 ==========
    try
        long_name  = ncreadatt(ncpath,varname,'long_name');
    catch
        long_name  = '无';
    end
    try
        comment    = ncreadatt(ncpath,varname,'comment');
    catch
        comment    = '无';
    end
    try
        units      = ncreadatt(ncpath,varname,'units');
    catch
        units      = '无';
    end
    try
        missing_val= ncreadatt(ncpath,varname,'missing_value');
    catch
        missing_val= NaN;
    end
    try
        fill_val   = ncreadatt(ncpath,varname,'_FillValue');
    catch
        fill_val   = NaN;
    end
    % 存入属性表
    attr_all(i).file      = fname;
    attr_all(i).model     = model;
    attr_all(i).table_id  = table_id;
    attr_all(i).long_name = long_name;
    attr_all(i).comment   = comment;
    attr_all(i).units     = units;
    attr_all(i).missing_value = missing_val;
    attr_all(i).FillValue = fill_val;
    
    % ========== 4.3 读取数据并处理缺失值 ==========
    cSoil = ncread(ncpath,varname);
    lon   = ncread(ncpath,'lon');
    lat   = ncread(ncpath,'lat');
    time  = ncread(ncpath,'time');
    
    hl = size(cSoil);
    attr_all(i).size = hl; 
    attr_all(i).lon = lon; 
    attr_all(i).lat = lat; 
    
    % 缺失值替换为NaN
    if ~isnan(fill_val)
        cSoil(cSoil==fill_val) = NaN;
    end
    if ~isnan(missing_val)
        cSoil(cSoil==missing_val) = NaN;
    end
    cSoil(cSoil==0) = NaN;
    
    % ========== 4.5 Emon月尺度→年平均；Eyr直接使用（支持多文件拼接）==========
    if strcmp(table_id,'Emon') || strcmp(table_id,'Lmon')
        % Emon：每月数据，每年12个月平均
        years = t_start:t_end;
        ny = length(years);
        cSoil_year = zeros(hl(1),hl(2),ny);
    
        % 计算起始累计月份（适配你定义的 m_start）
        total_month_before = 0;
        for y = 1:ny
            % 确定当前年份的起始月、结束月
            if y == 1
                s = m_start;
                e = 12;
            elseif y == ny
                s = 1;
                e = n_end;
            else
                s = 1;
                e = 12;
            end

            n_month = e - s + 1;
            idx = total_month_before + (1:n_month);
            idx(idx > size(cSoil, 3)) = [];

            cSoil_year(:,:,y) = mean(cSoil(:,:,idx),3,'omitnan');
        end

        % ============== 修复：同一个模式多个文件自动拼接 ==============
        if isfield(cSoil_year_all, model)
            % 模式已存在 → 拼接年份 & 数据
            cSoil_year_all.(model).data = cat(3, cSoil_year_all.(model).data, cSoil_year);
            cSoil_year_all.(model).years = [cSoil_year_all.(model).years, years];
        else
            % 模式不存在 → 初始化
            cSoil_year_all.(model).data = cSoil_year;
            cSoil_year_all.(model).years = years;
        end

    elseif strcmp(table_id,'Eyr')
        % Eyr：本身是年数据，直接使用
        years = t_start:t_end;
        cSoil_year = cSoil;

        % ============== 修复：同一个模式多个文件自动拼接 ==============
        if isfield(cSoil_year_all, model)
            cSoil_year_all.(model).data = cat(3, cSoil_year_all.(model).data, cSoil_year);
            cSoil_year_all.(model).years = [cSoil_year_all.(model).years, years];
        else
        cSoil_year_all.(model).data = cSoil_year;
        cSoil_year_all.(model).years = years;
        end
    end

    % 这里也要更新为拼接后的完整年份
    attr_all(i).years = years;
end

% 5. 模式统计 & 属性输出
fprintf('\n===== 模式统计结果 =====\n');
fprintf('共 %d 个气候模式：\n',length(model_list));
disp(model_list');

% 导出属性到命令行+表格
attr_table = struct2table(attr_all);
disp('===== 各文件 cSoil 属性汇总 =====');
disp(attr_table);
% 可选：导出
% save('H:\数据库\CMIP6_SSP585_Soil_Carbon_data\cSoil_属性汇总.mat', 'attr_table', '-v6');
% save('H:\数据库\CMIP6_SSP585_Soil_Carbon_data\cSoilSlow_属性汇总.mat', 'attr_table', '-v6');

% 6. 多模式集合平均（重采样到最高分辨率 + 格点级同年平均）
% -------------------------- 6.1 找到最高分辨率网格（目标网格） --------------------------
% 遍历所有文件，找到 512×256 的最高分辨率经纬度（你的截图里最大就是512×256）
max_nlon = 0; max_nlat = 0;
lon_target = []; lat_target = [];
for i = 1:length(attr_all)
    sz = attr_all(i).size;
    nlon_i = sz(1); nlat_i = sz(2);
    if nlon_i > max_nlon && nlat_i > max_nlat
        max_nlon = nlon_i;
        max_nlat = nlat_i;
        lon_target = attr_all(i).lon;
        lat_target = attr_all(i).lat;
    end
end
fprintf('最高分辨率：%d×%d，作为重采样目标网格\n',max_nlon,max_nlat);

% -------------------------- 6.2 提取所有年份 --------------------------
all_years = [];
for m = 1:length(model_list)
    mod = model_list{m};
        % 安全获取该模式的年份（自动转列向量）
    idx_mod = strcmp({attr_all.model}, mod);
    years_mod = attr_all(idx_mod).years;
    years_mod = years_mod(:);  % 强制转列向量，解决维度报错！
    
    all_years = [all_years; years_mod];
end
unique_years = unique(all_years);
unique_years = sort(unique_years);
ny_total = length(unique_years);
nlon = length(lon_target);
nlat = length(lat_target);

% 初始化集合平均结果：lon×lat×year
cSoil_multi_mean = zeros(nlon, nlat, ny_total);
cSoil_multi_sd = zeros(nlon, nlat, ny_total);
cSoil_multi_p05 = zeros(nlon, nlat, ny_total);
cSoil_multi_p95 = zeros(nlon, nlat, ny_total);

% 存储每个模式的距平时间序列
cSoil_model_anom = {};   % 每个模式的逐年数据 [year, cSoil_mean]
% -------------------------- 6.3 逐年重采样 + 多模式格点平均 --------------------------
for yi = 1:ny_total
    y = unique_years(yi);
    % 存储当前年份所有模式重采样后的数据
    mod_data_list = [];
    
    for m = 1:length(model_list)
        mod = model_list{m};
        
        % 从拼接好的结构体读取数据和年份
        data_mod = cSoil_year_all.(mod).data;
        years_mod = cSoil_year_all.(mod).years;
        
        % ===================== 关键：找到该模式 2015 年的基准值 =====================
        idx_2015 = find(years_mod == 2015);
        if isempty(idx_2015)
            continue; % 无2015年，跳过
        end
        data_2015 = data_mod(:,:,idx_2015); % 2015基准场
        % ===================== 新加的 =====================
        
        % 查找年份索引
        idx_y = find(years_mod == y);
        if isempty(idx_y)
            continue;
        end

        data_y = data_mod(:,:,idx_y); % lon×lat
        
        % ===================== 计算距平：当前年份 - 2015年 =====================
        data_y = data_y - data_2015;
        % ===================== 新加的 =====================
        
        % 读取该模式原始经纬度
        idx_mod = strcmp({attr_all.model}, mod);
        lon_mod = attr_all(idx_mod).lon;
        lat_mod = attr_all(idx_mod).lat;
        
        % ========== 关键：重采样到最高分辨率网格（双线性插值） ==========
        % 生成原始二维网格
        [X, Y] = meshgrid(lon_mod, lat_mod);
        % 生成目标二维网格
        [Xq, Yq] = meshgrid(lon_target, lat_target);

        % 插值（data_y 是 lon×lat）
        data_resamp = interp2(X, Y, data_y', Xq, Yq, 'linear', NaN);
        
        % ===================== ? 存储每个模式的距平数据 =====================
        if ~isfield(cSoil_model_anom, mod)
            % 第一次：初始化
            cSoil_model_anom.(mod).data = data_resamp;
            cSoil_model_anom.(mod).years = y;
        else
            % 已有数据：按年份拼接
            cSoil_model_anom.(mod).data = cat(3, cSoil_model_anom.(mod).data, data_resamp);
            cSoil_model_anom.(mod).years = [cSoil_model_anom.(mod).years, y];
        end
        
        mod_data_list = cat(3, mod_data_list, data_resamp);
    end
    
    % 同一年份、同一格点：所有模式平均
    cSoil_multi_mean(:,:,yi) = mean(mod_data_list, 3, 'omitnan')';
    cSoil_multi_sd(:,:,yi) = std(mod_data_list, 0, 3, 'omitnan')';
    cSoil_multi_p05(:,:,yi)  = prctile(mod_data_list, 5, 3)';   % 5%分位数
    cSoil_multi_p95(:,:,yi)  = prctile(mod_data_list, 95, 3)';  % 95%分位数
end

save('H:\数据库\CMIP6_SSP585_Soil_Carbon_data\碳库距平后平均数据\cSoilMedium_多模式结果.mat', ...
    'cSoil_multi_mean','cSoil_multi_sd','cSoil_multi_p05','cSoil_multi_p95',...
    'cSoil_model_anom', ...  % 每个模式的距平场
    'lon_target','lat_target','unique_years','model_list',...
    'max_nlat', 'max_nlon',...
    'attr_table', '-v6');

figure(1)
image(cSoil_multi_mean(:,:,1))

% 7. 输出结果、绘图、保存
% 1. 输出全球逐年平均一维时间序列
global_ts = squeeze(mean(mean(cSoil_multi_mean,1,'omitnan'),2,'omitnan'));
result = [unique_years, global_ts];
disp('===== 多模式集合平均全球年平均时间序列 [年份, cSoil全球年平均] =====');
disp(result);

% 2. 绘制全球平均时间序列
figure('Color','w');
plot(unique_years, global_ts, '-o','LineWidth',1.5,'MarkerSize',4);
xlabel('年份');
ylabel('cSoil 全球年平均土壤碳库');
title('CMIP6 SSP585 多模式集合平均土壤碳库时间序列');
grid on;

% 3. 绘制第一年空间分布图示例
figure('Color','w');
pcolor(lon_target, lat_target, cSoil_multi_mean(:,:,15)');
shading interp; colorbar;
xlabel('经度'); ylabel('纬度');
title('多模式集合平均 第一年土壤碳库空间分布');

result = [unique_years, global_ts];
disp('===== 多模式集合平均全球年平均时间序列 [年份, cSoil全球年平均] =====');
disp(result);

% 绘图
figure('Color','w');

% ---------------- 绘制 5-95% 灰色阴影 ----------------
fill([unique_years; flip(unique_years)], ...   % 闭合X轴
     [global_p05; flip(global_p95)], ...       % 闭合Y轴
     [0.9 0.9 0.9], ...                        % 浅灰色
     'EdgeColor', 'none', 'FaceAlpha', 0.5);   % 透明阴影
hold on;

% ---------------- 绘制集合平均曲线 ----------------
plot(unique_years, global_ts, '-o','LineWidth',1.5,'MarkerSize',4, 'Color', [0.1 0.4 0.8]);

xlabel('年份');
ylabel('cSoil 全球年平均土壤碳库');
title('CMIP6 SSP585 多模式集合平均土壤碳库时间序列');
grid on;
legend('5–95% 范围','集合平均','Location','best');


% 8. 将多模式集合平均保存为 带地理信息的 TIFF 文件 (ArcGIS直接用)
type='Slow';
load(fullfile('H:\数据库\CMIP6_SSP585_Soil_Carbon_data\碳库距平后平均数据', ['cSoil',type,'_多模式结果.mat']))

% 地理参考信息设置 (CMIP6 标准: 全球等经纬度投影)
R = georasterref(...
    'RasterSize', [max_nlat, max_nlon], ...
    'Latlim', [double(min(lat_target)), double(max(lat_target))], ...
    'Lonlim', [double(min(lon_target)), double(max(lon_target))], ...
    'ColumnsStartFrom', 'north');
yy=[2026,2039,2060,2078,2094];
% 计算相对于2015年的变化
% cSoil_multi_mean = cSoil_multi_mean-cSoil_multi_mean(:,:,1);
for ii=1:5
    tif_name = fullfile('H:\数据库\CMIP6_SSP585_Soil_Carbon_data\TIFF_输出2', ['cSoil',type,sprintf('_mean_%d_change.tif', yy(ii))]);
    geotiffwrite(tif_name, flipud(cSoil_multi_mean(:,:,yy(ii)-2015+1)'), R, 'CoordRefSysCode', 4326); 
end

%% 计算碳源汇窗口的空间分布
% 替换为相对重要性和敏感性的加权
% 使用更新的碳源汇概率曲线数据
load('H:\数据库\碳源汇预测曲线3\更新了碳源汇概率曲线.mat')

ssp=1; %1.5
% ssp=2; %2.0
% ssp=3; %3.0
% ssp=4; %4.0
% ssp=5; %5.0
env_res = [env_results{1}.ssp_results];
bin_edges = env_res(ssp).bin_edges;
mat_g2_lower = bin_edges(1:50)';
mat_g2_upper = bin_edges(2:51)';
mat_g2_col5 = env_res(ssp).prob_sink';
mat_g3_lower = bin_edges(1:50)';
mat_g3_upper = bin_edges(2:51)';
mat_g3_col5 = env_res(ssp).prob_source';
env_res = [env_results{2}.ssp_results];
bin_edges = env_res(ssp).bin_edges;
map_g2_lower = bin_edges(1:50)';
map_g2_upper = bin_edges(2:51)';
map_g2_col5 = env_res(ssp).prob_sink';
map_g3_lower = bin_edges(1:50)';
map_g3_upper = bin_edges(2:51)';
map_g3_col5 = env_res(ssp).prob_source';
env_res = [env_results{3}.ssp_results];
bin_edges = env_res(ssp).bin_edges;
et_g2_lower = bin_edges(1:50)';
et_g2_upper = bin_edges(2:51)';
et_g2_col5 = env_res(ssp).prob_sink';
et_g3_lower = bin_edges(1:50)';
et_g3_upper = bin_edges(2:51)';
et_g3_col5 = env_res(ssp).prob_source';
env_res = [env_results{4}.ssp_results];
bin_edges = env_res(ssp).bin_edges;
ai_g2_lower = bin_edges(1:50)'*10000;
ai_g2_upper = bin_edges(2:51)'*10000;
ai_g2_col5 = env_res(ssp).prob_sink';
ai_g3_lower = bin_edges(1:50)'*10000;
ai_g3_upper = bin_edges(2:51)'*10000;
ai_g3_col5 = env_res(ssp).prob_source';
env_res = [env_results{5}.ssp_results];
bin_edges = env_res(ssp).bin_edges;
npp_g2_lower = bin_edges(1:50)';
npp_g2_upper = bin_edges(2:51)';
npp_g2_col5 = env_res(ssp).prob_sink';
npp_g3_lower = bin_edges(1:50)';
npp_g3_upper = bin_edges(2:51)';
npp_g3_col5 = env_res(ssp).prob_source';

load('H:\数据库\碳源汇预测曲线3\更新了碳源汇概率曲线2_POCMAOCTOC.mat')
env_res = [env_results{6}.ssp_results];
bin_edges = env_res(ssp).bin_edges;
POC_g2_lower = bin_edges(1:50)';
POC_g2_upper = bin_edges(2:51)';
POC_g2_col5 = env_res(ssp).prob_sink';
POC_g3_lower = bin_edges(1:50)';
POC_g3_upper = bin_edges(2:51)';
POC_g3_col5 = env_res(ssp).prob_source';

env_res = [env_results{7}.ssp_results];
bin_edges = env_res(ssp).bin_edges;
MAOC_g2_lower = bin_edges(1:50)';
MAOC_g2_upper = bin_edges(2:51)';
MAOC_g2_col5 = env_res(ssp).prob_sink';
MAOC_g3_lower = bin_edges(1:50)';
MAOC_g3_upper = bin_edges(2:51)';
MAOC_g3_col5 = env_res(ssp).prob_source';

% 全球环境值
load('H:\数据库\随机森林相对重要性\全球所有点的历史环境值.mat')

indir='H:\数据库\随机森林相对重要性\未来情景的预测值\';
load([indir, 'SSP126_2140的环境值.mat']);
info = geotiffinfo('F:\参考蒸散量ET0数据\gaez-v5-data\temp.tif');
[temp2, R] = geotiffread('F:\参考蒸散量ET0数据\gaez-v5-data\temp.tif');

param_SSP_tif30 = single(zeros([size(temp2),7])); 
param_SSP_tif30(:) = NaN;
for ii=1:length(envi_SSP126_2140)
    ii
    if any(isnan(new_ENVI_hist(ii,:))) 
        continue
    end
    
    % 构造比较矩阵（m×n）
    comp_matrix = (new_ENVI_hist(ii,3) > mat_g2_lower) & (new_ENVI_hist(ii,3) <= mat_g2_upper);
    % 找到第一个匹配的行索引
    if isempty(find(comp_matrix, 1))
        continue
    end
    mat_g2 = mat_g2_col5(comp_matrix);
    
    comp_matrix = (new_ENVI_hist(ii,3) > mat_g3_lower) & (new_ENVI_hist(ii,3) <= mat_g3_upper);
    % 找到第一个匹配的行索引
    if isempty(find(comp_matrix, 1))
        continue
    end
    mat_g3 = mat_g3_col5(comp_matrix);
    
    comp_matrix = (new_ENVI_hist(ii,4) > map_g2_lower) & (new_ENVI_hist(ii,4) <= map_g2_upper);
    % 找到第一个匹配的行索引
    if isempty(find(comp_matrix, 1))
        continue
    end
    map_g2 = map_g2_col5(comp_matrix);
    
    comp_matrix = (new_ENVI_hist(ii,4) > map_g3_lower) & (new_ENVI_hist(ii,4) <= map_g3_upper);
    % 找到第一个匹配的行索引
    if isempty(find(comp_matrix, 1))
        continue
    end
    map_g3 = map_g3_col5(comp_matrix);
    
    comp_matrix = (new_ENVI_hist(ii,5) > et_g2_lower) & (new_ENVI_hist(ii,5) <= et_g2_upper);
    % 找到第一个匹配的行索引
    if isempty(find(comp_matrix, 1))
        continue
    end
    et_g2 = et_g2_col5(comp_matrix);
    
    comp_matrix = (new_ENVI_hist(ii,5) > et_g3_lower) & (new_ENVI_hist(ii,5) <= et_g3_upper);
    % 找到第一个匹配的行索引
    if isempty(find(comp_matrix, 1))
        continue
    end
    et_g3 = et_g3_col5(comp_matrix);
    
    comp_matrix = (new_ENVI_hist(ii,6) > ai_g2_lower) & (new_ENVI_hist(ii,6) <= ai_g2_upper);
    % 找到第一个匹配的行索引
    if isempty(find(comp_matrix, 1))
        continue
    end
    ai_g2 = ai_g2_col5(comp_matrix);
    
	comp_matrix = (new_ENVI_hist(ii,6) > ai_g3_lower) & (new_ENVI_hist(ii,6) <= ai_g3_upper);
    % 找到第一个匹配的行索引
    if isempty(find(comp_matrix, 1))
        continue
    end
    ai_g3 = ai_g3_col5(comp_matrix);
    
    comp_matrix = (new_ENVI_hist(ii,7) > npp_g2_lower) & (new_ENVI_hist(ii,7) <= npp_g2_upper);
    % 找到第一个匹配的行索引
    if isempty(find(comp_matrix, 1))
        continue
    end
    npp_g2 = npp_g2_col5(comp_matrix);

	comp_matrix = (new_ENVI_hist(ii,7) > npp_g3_lower) & (new_ENVI_hist(ii,7) <= npp_g3_upper);
    % 找到第一个匹配的行索引
    if isempty(find(comp_matrix, 1))
        continue
    end
    npp_g3 = npp_g3_col5(comp_matrix);
    
    g2 = [mat_g2,map_g2,et_g2,ai_g2,npp_g2] - [mat_g3,map_g3,et_g3,ai_g3,npp_g3];
    param_SSP_tif30(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,1:5) = g2;
end

indir_npp = 'H:\数据库\NPP_NDVI_MODIS\Global_Mean_NPP_2003_2023_MODIS重采样0.0833.tif';
% 1.2 读取NPP
[meanNPP, R_npp] = geotiffread(indir_npp);
% 基准栅格行列数
rows_npp = size(meanNPP,1);
cols_npp = size(meanNPP,2);
disp(['基准栅格尺寸：',num2str(rows_npp),' × ',num2str(cols_npp)]);
% 从你R_npp显示的属性提取边界、步长
lon_min = R_npp.LongitudeLimits(1);
lon_max = R_npp.LongitudeLimits(2);
dx = R_npp.CellExtentInLongitude;

lat_min = R_npp.LatitudeLimits(1);
lat_max = R_npp.LatitudeLimits(2);
dy = R_npp.CellExtentInLatitude;

% 构建目标经纬度网格
lon_vec = lon_min : dx : (lon_max - dx);
lat_vec = lat_max : -dy : (lat_min + dy);
[lon_target, lat_target] = meshgrid(lon_vec, lat_vec);

% MAOC  这个没变  是使用重采样后的BD 单位转换后的结果
[MAOC0_30, R_maoc] = geotiffread('H:\数据库\MAOC0-30cm.tif');  

% POC 这个变了 得将20cm间隔的做平均 
[POC, R_poc] = geotiffread('H:\数据库\估计的POC0-100cm每隔20cm一层.tif');   % 0-100cm 以20cm间隔
POC0_30 = mean(cat(3,POC(:,:,1),POC(:,:,1),POC(:,:,2)), 3,'omitnan');
clear POC
MAOC0_30 = double(MAOC0_30);
POC0_30 = double(POC0_30);

% 构建MAOC原始坐标
lon_maoc_min = R_maoc.LongitudeLimits(1);
lon_maoc_max = R_maoc.LongitudeLimits(2);
dx_maoc = R_maoc.CellExtentInLongitude;
lat_maoc_min = R_maoc.LatitudeLimits(1);
lat_maoc_max = R_maoc.LatitudeLimits(2);
dy_maoc = R_maoc.CellExtentInLatitude;

lon_maoc_vec = lon_maoc_min : dx_maoc : (lon_maoc_max - dx_maoc);
lat_maoc_vec = lat_maoc_max : -dy_maoc : (lat_maoc_min + dy_maoc);
[lon_maoc, lat_maoc] = meshgrid(lon_maoc_vec, lat_maoc_vec);

% 构建POC原始坐标
lon_poc_min = R_poc.LongitudeLimits(1);
lon_poc_max = R_poc.LongitudeLimits(2);
dx_poc = R_poc.CellExtentInLongitude;
lat_poc_min = R_poc.LatitudeLimits(1);
lat_poc_max = R_poc.LatitudeLimits(2);
dy_poc = R_poc.CellExtentInLatitude;

lon_poc_vec = lon_poc_min : dx_poc : (lon_poc_max - dx_poc);
lat_poc_vec = lat_poc_max : -dy_poc : (lat_poc_min + dy_poc);
[lon_poc, lat_poc] = meshgrid(lon_poc_vec, lat_poc_vec);

% 展平
lon_maoc_flat = lon_maoc(:);
lat_maoc_flat = lat_maoc(:);
maoc_flat = MAOC0_30(:);

% 用scatteredInterpolant，大数据比griddata更快，域外=NaN
F_maoc = scatteredInterpolant(lon_maoc_flat, lat_maoc_flat, maoc_flat, 'linear','none');
MAOC_resamp = F_maoc(lon_target, lat_target);

lon_poc_flat = lon_poc(:);
lat_poc_flat = lat_poc(:);
poc_flat = POC0_30(:);

F_poc = scatteredInterpolant(lon_poc_flat, lat_poc_flat, poc_flat, 'linear','none');
POC_resamp = F_poc(lon_target, lat_target);

for ii=1:2160
    for jj=1:4320
        ii
        if isnan(POC_resamp(ii,jj)) || isnan(MAOC_resamp(ii,jj))
            continue
        end
        comp_matrix = (POC_resamp(ii,jj) > POC_g2_lower) & (POC_resamp(ii,jj) <= POC_g2_upper);
        % 找到第一个匹配的行索引
        if isempty(find(comp_matrix, 1))
            continue
        end
        POC_g2 = POC_g2_col5(comp_matrix);

        comp_matrix = (POC_resamp(ii,jj) > POC_g3_lower) & (POC_resamp(ii,jj) <= POC_g3_upper);
        % 找到第一个匹配的行索引
        if isempty(find(comp_matrix, 1))
            continue
        end
        POC_g3 = POC_g3_col5(comp_matrix);
        
        param_SSP_tif30(ii,jj,6) = POC_g2-POC_g3;
    
        comp_matrix = (MAOC_resamp(ii,jj) > MAOC_g2_lower) & (MAOC_resamp(ii,jj) <= MAOC_g2_upper);
        % 找到第一个匹配的行索引
        if isempty(find(comp_matrix, 1))
            continue
        end
        MAOC_g2 = MAOC_g2_col5(comp_matrix);

        comp_matrix = (MAOC_resamp(ii,jj) > MAOC_g3_lower) & (MAOC_resamp(ii,jj) <= MAOC_g3_upper);
        % 找到第一个匹配的行索引
        if isempty(find(comp_matrix, 1))
            continue
        end
        MAOC_g3 = MAOC_g3_col5(comp_matrix);
        
        param_SSP_tif30(ii,jj,7) = MAOC_g2-MAOC_g3;
        
    end
end



% 权重是通过 H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\RF结果6\相对重要性.xlsx 中的数据求平均计算的
aa=(param_SSP_tif30(:,:,1)*0.45+param_SSP_tif30(:,:,2)*0.05+param_SSP_tif30(:,:,3)*0.03+ ...
    param_SSP_tif30(:,:,4)*0.04+param_SSP_tif30(:,:,5)*0.5);

% geotiffwrite('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\碳源汇预测曲线3\全球30_5.0TOC碳汇窗口指数new加权.tif',aa,R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag)
% geotiffwrite('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\碳源汇预测曲线3\全球30_4.0TOC碳汇窗口指数new加权.tif',aa,R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag)
% geotiffwrite('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\碳源汇预测曲线3\全球30_3.0TOC碳汇窗口指数new加权.tif',aa,R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag)
% geotiffwrite('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\碳源汇预测曲线3\全球30_2.0TOC碳汇窗口指数new加权.tif',aa,R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag)
geotiffwrite('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\碳源汇预测曲线3\全球30_1.5TOC碳汇窗口指数new加权.tif',aa,R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag)

% 用斜率计算敏感性
[aa15, ~] = geotiffread('H:\数据库\碳源汇预测曲线3\全球30_1.5TOC碳汇窗口指数new加权.tif'); 
[aa20, ~] = geotiffread('H:\数据库\碳源汇预测曲线3\全球30_2.0TOC碳汇窗口指数new加权.tif'); 
[aa30, ~] = geotiffread('H:\数据库\碳源汇预测曲线3\全球30_3.0TOC碳汇窗口指数new加权.tif'); 
[aa40, ~] = geotiffread('H:\数据库\碳源汇预测曲线3\全球30_4.0TOC碳汇窗口指数new加权.tif'); 
[aa50, ~] = geotiffread('H:\数据库\碳源汇预测曲线3\全球30_5.0TOC碳汇窗口指数new加权.tif'); 

% 1. 自变量设置
x = [1.5, 2.0, 3.0, 4.0, 5.0];
n = length(x);

% 2. 拼接三维数组 (行,列,层数)
Y = cat(3, aa15, aa20, aa30, aa40, aa50);

% 3. 预分配斜率结果矩阵
[rowNum, colNum, ~] = size(Y);
slopeMap = nan(rowNum, colNum, 'single');

% 4. 逐像元循环计算斜率
for ii = 1:rowNum
    for jj = 1:colNum
        y_vec = squeeze(Y(ii,jj,:));
        
        % 存在NaN则该位置保持NaN
        if any(isnan(y_vec))
            continue;
        end
        
        % 一阶线性拟合 p(1)=斜率k, p(2)=截距b
        p = polyfit(x, y_vec', 1);
%         k = p(1);
%         b = p(2);
        slopeMap(ii,jj) = p(1);
    end
end

geotiffwrite('H:\数据库\碳源汇预测曲线3\全球30TOC碳汇窗口指数new加权线性拟合敏感性.tif',slopeMap,R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag)

%% 计算碳源汇概率曲线，新修改的。
clear; clc; close all;
% 路径配置
indir_env = 'F:\黄土高原\';
indir_ssp = 'H:\数据库\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\';     
% ====================== 加载所有环境变量数据 ======================
% 1.1 加载ET0/AI/MAT/MAP
load([indir_env, '参考蒸散量ET0数据\gaez-v5-data\历史数据\1980-2022年球ET0平均值.mat']); % meanET0
load([indir_env, '干旱指数\用MAP计算的AI\1980-2022年球AI平均值.mat']); % meanAI1
load([indir_env, 'ERA5气候数据\20251111下载的新的\1980-2024年球年平均气温平均值.mat']); % meanMAT
load([indir_env, 'ERA5气候数据\20251111下载的新的\1980-2024年球年降水量平均值.mat']); % meanMAP
indir_npp = 'H:\数据库\NPP_NDVI_MODIS\Global_Mean_NPP_2003_2023_MODIS重采样0.0833.tif';
% 1.2 读取NPP
[meanNPP, R_npp] = geotiffread(indir_npp);

% 1.3 统一预处理MAT/MAP（维度匹配）
% MAT预处理
result=meanMAT;
result = result(:,1:1800,:);
result = rot90(result, -1); 
result = flip(result, 2);
result= [result(:,1801:3600,:),result(:,1:1800,:)]; 
result2 = imresize(result, [2160 4320], 'nearest');
meanMAT=result2; clear result result2;

% MAP预处理
result=meanMAP;
result = result(:,1:1800,:);
result = rot90(result, -1); 
result = flip(result, 2);
result= [result(:,1801:3600,:),result(:,1:1800,:)]; 
result2 = imresize(result, [2160 4320], 'nearest');
meanMAP=result2; clear result result2;

% ====================== 定义核心参数 ======================
% SSP情景配置
ssp_list = {'1.5', '2.0', '3.0', '4.0', '5.0'};
ssp_files = {
    [indir_ssp, '全球碳库SSP预测30_1.5C增量TOC.tif'], ...
    [indir_ssp, '全球碳库SSP预测30_2.0C增量TOC.tif'], ...
    [indir_ssp, '全球碳库SSP预测30_3.0C增量TOC.tif'], ...
    [indir_ssp, '全球碳库SSP预测30_4.0C增量TOC.tif'], ...
    [indir_ssp, '全球碳库SSP预测30_5.0C增量TOC.tif'] ...
};
% 线型：5 种
ssp_line_style = {'-', '-', ':', '-.', '--'}; 

% 碳源：蓝色系（冷色）
ssp_color_source = [
    0.0 0.2 0.6;
    0.1 0.3 0.8;
    0.2 0.4 1.0;
    0.4 0.6 1.0;
    0.6 0.8 1.0;
];  

% 碳汇：橙黄色系（暖色）
ssp_color_sink = [
    0.6 0.3 0.0;
    0.8 0.4 0.0;
    1.0 0.5 0.0;
    1.0 0.7 0.2;
    1.0 0.8 0.4;
];   
env_vars = { ...
    % 变量内部名, 显示名, 数据, 过滤最小值, 过滤最大值
    {'MAT', '年平均气温(℃)', meanMAT, -14, 29}, ...
    {'MAP', '年降水量(mm)', meanMAP, 200, 3000}, ...
    {'ET0', '参考蒸散量(mm)', meanET0, 160, 2800}, ...
    {'AI', '干旱指数', meanAI1/10000, 600/10000, 29000/10000}, ...
    {'NPP', '净初级生产力(kg m-2 a-1)', meanNPP, 0.006, 1.57}, ...
};

% 计算参数
n_bins = 50; % 分箱数
window_width = 7; % 平滑窗宽

% 3. 批量计算所有环境变量×SSP情景的概率
env_results = cell(length(env_vars), 1); % 存储每个环境变量的结果

for env_idx = 1:length(env_vars)
    fprintf('\n=====================================\n');
    fprintf('开始计算：%s\n', env_vars{env_idx}{2}); % 取显示名
    fprintf('=====================================\n');
    
    % 获取当前环境变量信息
    env_name = env_vars{env_idx}{1};
    env_display = env_vars{env_idx}{2};
    env_data = env_vars{env_idx}{3};
    env_min = env_vars{env_idx}{4};
    env_max = env_vars{env_idx}{5};
    
    % 初始化当前环境变量的SSP结果（元胞数组）
    env_ssp_results = cell(length(ssp_list), 1);
    
    % 遍历每个SSP情景
    for ssp_idx = 1:length(ssp_list)
        fprintf('\n----- 计算 %s 情景 -----\n', ssp_list{ssp_idx});
        % 读取SSP的TOC数据
        [TOC, ~] = geotiffread(ssp_files{ssp_idx});
        % 计算概率
        [bin_centers, prob_source, prob_sink, var_min, var_max, bin_edges] = calc_env_prob2(...
            env_data, meanMAT, meanMAP, meanAI1/10000, meanET0, meanNPP, TOC, n_bins);
        % 移动平均平滑
        prob_source_smooth = smoothdata(prob_source, 'movmean', window_width);
        prob_sink_smooth = smoothdata(prob_sink, 'movmean', window_width);
        % 存储结果（结构体赋值到元胞数组）
        env_ssp_results{ssp_idx} = struct( ...
            'name', ssp_list{ssp_idx}, ...
            'bin_edges', bin_edges, ...
            'bin_centers', bin_centers, ...
            'prob_source', prob_source_smooth, ...
            'prob_sink', prob_sink_smooth, ...
            'var_min', var_min, ...
            'var_max', var_max ...
        );
    
    end
    
    % 存储当前环境变量的所有SSP结果（结构体）
    env_results{env_idx} = struct( ...
        'name', env_name, ...
        'display_name', env_display, ...
        'ssp_results', env_ssp_results ... % 元胞数组字段
    );
end
save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\碳源汇预测曲线3\更新了碳源汇概率曲线.mat', 'env_results', '-v6');

% 1. 创建总画布（统一设置大小和背景）
figure('Color','w','Position',[100, 100, 1800, 1200]); % 增大画布适配2行3列

% 2. 循环绘制每个子图（2行3列，共6个子图）
for env_idx = 1:length(env_results)
    % 定位到当前子图（2行3列，第env_idx个位置）
    subplot(3, 3, env_idx);
    hold on; % 保持当前子图，用于叠加多条曲线
    
    % 获取当前环境变量信息
    env_res = env_results{env_idx};
    env_name = env_res.name;
    ssp_res_list = [env_res.ssp_results]; % 转为结构体数组
    
    % 绘制碳源曲线（红色系）
    for ssp_idx = 1:length(ssp_res_list)
        ssp_res = ssp_res_list(ssp_idx);
        plot(ssp_res.bin_centers, ssp_res.prob_source, ...
            'Color', ssp_color_source(ssp_idx,:), ...
            'LineStyle', ssp_line_style{ssp_idx}, ...
            'LineWidth', 2, ...
            'DisplayName', [ssp_res.name, ' - Carbon source']);
    end
    
    % 绘制碳汇曲线（绿色系）
    for ssp_idx = 1:length(ssp_res_list)
        ssp_res = ssp_res_list(ssp_idx);
        plot(ssp_res.bin_centers, ssp_res.prob_sink, ...
            'Color', ssp_color_sink(ssp_idx,:), ...
            'LineStyle', ssp_line_style{ssp_idx}, ...
            'LineWidth', 2, ...
            'DisplayName', [ssp_res.name, ' - Carbon sink']);
    end
    
    % 子图美化
    xlabel(env_name, 'FontSize',10,'FontWeight','bold');
    ylabel('加权概率（0-1）','FontSize',10,'FontWeight','bold');
%     title(['SSP情景下', env_display, '碳源/碳汇加权概率'], ...
%         'FontSize',11,'FontWeight','bold');
    grid on; grid minor;
    box on;
    ylim([0, 1.05]);
    set(gca, 'GridAlpha',0.3, 'MinorGridAlpha',0.1);
    set(gca, 'FontSize',9, 'FontName','SimHei'); % 支持中文
    
    % 仅在第一个子图显示图例（避免重复），或统一放在画布底部
    if env_idx == 1
        legend('Location','best','FontSize',8);
    end
end


% ====================== 新增：初始化数据存储（改用元胞数组+统一维度） ======================
curve_data_cell = {};  % 用元胞数组暂存，最后统一转换为table
row_idx = 1;

for env_idx = 1:length(env_results)
    % 获取当前环境变量信息
    env_res = env_results{env_idx};
    env_name = env_res.name;
    ssp_res_list = [env_res.ssp_results]; % 转为结构体数组
    
    % 遍历SSP情景，提取碳源/碳汇数据
    for ssp_idx = 1:length(ssp_res_list)
        ssp_res = ssp_res_list(ssp_idx);
        ssp_name = ssp_res.name;
        
        % -------- 统一获取数据长度（避免维度不一致） --------
        n_points = length(ssp_res.bin_centers);
        
        % -------- 提取碳源数据（强制转为字符数组，统一长度） --------
        env_name_arr = repmat(char(env_name), n_points, 1);       % 转为字符数组
        ssp_name_arr = repmat(char(ssp_name), n_points, 1);       % 转为字符数组
        curve_type_source = repmat(char('Carbon source'), n_points, 1);  % 统一字符长度
        
        % 构建碳源数据行（元胞数组存储，避免直接拼接table）
        for i = 1:n_points
            curve_data_cell{row_idx, 1} = env_name_arr(i, :);    % 环境变量名
            curve_data_cell{row_idx, 2} = ssp_name_arr(i, :);    % SSP情景
            curve_data_cell{row_idx, 3} = curve_type_source(i, :);% 曲线类型
            curve_data_cell{row_idx, 4} = double(ssp_res.bin_centers(i)); % X轴：bin中心
            curve_data_cell{row_idx, 5} = double(ssp_res.prob_source(i)); % Y轴：加权概率
            row_idx = row_idx + 1;
        end
        
        % -------- 提取碳汇数据（同上，统一格式） --------
        curve_type_sink = repmat(char('Carbon sink'), n_points, 1);
        for i = 1:n_points
            curve_data_cell{row_idx, 1} = env_name_arr(i, :);
            curve_data_cell{row_idx, 2} = ssp_name_arr(i, :);
            curve_data_cell{row_idx, 3} = curve_type_sink(i, :);
            curve_data_cell{row_idx, 4} = double(ssp_res.bin_centers(i));
            curve_data_cell{row_idx, 5} = double(ssp_res.prob_sink(i));
            row_idx = row_idx + 1;
        end
    end
end

% ====================== 转换为table并保存（解决维度不匹配） ======================
% 将元胞数组转换为table（自动适配所有行的维度）
curve_data_all = cell2table(curve_data_cell, 'VariableNames', {'env_name', 'ssp_name', 'curve_type', 'x', 'y'});

% ====================== 保存数据到本地 ======================
% 方式2：保存为CSV文件
writetable(curve_data_all, 'H:\数据库\碳源汇预测曲线3\更新了碳源汇概率曲线POCMAOCTOC.csv', 'Encoding', 'UTF-8');
%% 提取典型区，不同过程和碳库随着增温的变化

% 读取参数值
indir='H:\数据库\随机森林相对重要性\未来情景的预测值\';
param_name = {'kbd_0_30_83', 'kl_0_30_83', 'kp_0_30_83', 'kq_0_30_83', 'CUE'};
yt = {'2140','4160','6180','8100'};
hh = 1837283;  %全变行数
param_SSP = single(zeros(hh,5,4)); % 只变一个变量预测值  9个变量 4个时间点
for ii =1:5
    for jj=1:4
        files = [indir, 'Global_',param_name{ii},'_SSP585_',yt{jj},'.mat'];
        temp = load(files);
        fieldNames = fieldnames(temp);
        originalVarName = fieldNames{1};
        temp2 = temp.(originalVarName);
        param_SSP(:,ii,jj)=temp2(:,3);
    end
end

indir='H:\数据库\随机森林相对重要性\未来情景的预测值\';
load([indir, 'SSP126_2140的环境值.mat']);
info = geotiffinfo('F:\参考蒸散量ET0数据\gaez-v5-data\temp.tif');
[lat,lon]=pix2latlon(info.RefMatrix,envi_SSP126_2140(:,1)+1,envi_SSP126_2140(:,2)+1); % 之前的索引从0开始的,在MATLAB中需要加1
latlon = [lat,lon];
clear lat lon

param_hist = single(zeros(hh,5));
% 读取数据
[kbd_0_30, R] = geotiffread('H:\数据库\0-30cm与30-100cm间隔分层0.008333度\kbd_0_30.tif');
[row, col] = geographicToDiscrete(R, latlon(:,1), latlon(:,2));  % 转换为栅格行列号
row(isnan(row)) = 1;
col(isnan(col)) = 1;
param_hist(:,1) = kbd_0_30(sub2ind(size(kbd_0_30), row, col)); %取出对应位置的值
clear kbd_0_30
[kl_0_30, ~] = geotiffread('H:\数据库\0-30cm与30-100cm间隔分层0.008333度\kl_0_30.tif');
param_hist(:,2) = kl_0_30(sub2ind(size(kl_0_30), row, col)); %取出对应位置的值
clear kl_0_30
[kp_0_30, ~] = geotiffread('H:\数据库\0-30cm与30-100cm间隔分层0.008333度\kp_0_30.tif');
param_hist(:,3) = kp_0_30(sub2ind(size(kp_0_30), row, col)); %取出对应位置的值
clear kp_0_30
[kq_0_30, ~] = geotiffread('H:\数据库\0-30cm与30-100cm间隔分层0.008333度\kq_0_30.tif');
param_hist(:,4) = kq_0_30(sub2ind(size(kq_0_30), row, col)); %取出对应位置的值
clear kq_0_30

% kbd 计算为kbd * MBC
% 只按0-30  30-100cm分层 单位 mgC kg-1Soil;
[mbc0_30, R2] = geotiffread( 'j:\数据备份-2023.12.18\global pattern\MBC30cm.tif');
mask2 = (mbc0_30 <= 0);
mbc0_30(mask2) = NaN;
mbc0_30 = mbc0_30/1000; % gC kg-1Soil;
clear mask2
[row, col] = geographicToDiscrete(R2, latlon(:,1), latlon(:,2));  % 转换为栅格行列号
row(isnan(row)) = 1;
col(isnan(col)) = 1;
param2_mbc0_30 = mbc0_30(sub2ind(size(mbc0_30), row, col)); %取出对应位置的值 
param_hist(:,1)=param_hist(:,1).*param2_mbc0_30;
clear mbc0_30 R2

% 加入CUE变量
% CUE 这个只用0-30cm单层的, 和MBC的分辨率一样
[CUE, R] = geotiffread('J:\数据备份-2023.12.18\global pattern\CUE_y.tif');
CUE(CUE<=0) = NaN;
[row, col] = geographicToDiscrete(R, latlon(:,1), latlon(:,2));  % 转换为栅格行列号
row(isnan(row)) = 1;
col(isnan(col)) = 1;
param2_cue = CUE(sub2ind(size(CUE), row, col)); %取出对应位置的值
param_hist(:,5) = param2_cue;
clear param2_cue CUE

kb = single(zeros([2160,4320,5]));
kl = single(zeros([2160,4320,5]));
kp = single(zeros([2160,4320,5]));
kq = single(zeros([2160,4320,5]));
cue = single(zeros([2160,4320,5]));
kb(:) = NaN;
kl(:) = NaN;
kp(:) = NaN;
kq(:) = NaN;
cue(:) = NaN;
for ii=1:length(envi_SSP126_2140)
    kb(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,1) = param_hist(ii,1);
    kl(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,1) = param_hist(ii,2);
    kp(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,1) = param_hist(ii,3);
    kq(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,1) = param_hist(ii,4);
    cue(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,1) = param_hist(ii,5);
    
    kb(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,2:5) = param_SSP(ii,1,:);
    kl(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,2:5) = param_SSP(ii,2,:);
    kp(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,2:5) = param_SSP(ii,3,:);
    kq(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,2:5) = param_SSP(ii,4,:);
    cue(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,2:5) = param_SSP(ii,5,:);
end

clearvars -except kb kl kp kq cue

% 典型区掩膜======================
[mask, ~] = geotiffread('H:\数据库\碳源汇预测曲线3\典型区\MASK2.tif');
 % 1 亚马孙  2 北方森林  3 北方森林北美  4 刚果  5 北极苔原 6 北极苔原北美
% 合并相同气候区
mask(mask==4)=1; % 1 亚马孙 4 刚果  ==热带雨林
mask(mask==3)=2; % 2 北方森林 3 北方森林北美  ==北方森林
mask(mask==5)=3; % 5 北极苔原 6 北极苔原北美  ==北极苔原
mask(mask==6)=3; % 5 北极苔原 6 北极苔原北美  ==北极苔原
 
% 1. 定义参数和区域列表
% 区域编号（按你说的顺序）
% regions = [1 2 3 4 5 6];
regions = [1 2 3]; % 合并相同气候区
% region_names = {'Amazon','Boreal_Eurasia','Boreal_NorthAmerica','Congo','Tundra_Eurasia','Tundra_NorthAmerica'};
region_names = {'Tropical_rainforest','Boreal_Eurasia','Tundra_Eurasia'}; % 合并相同气候区

% 变量名
carbon_vars = {'kb','kl','kp','kq','cue'};

% 时间点标签（对应第4维/第3维）
time_env = {'hist','2020-2040','2040-2060','2060-2080','2080-2100'};

% 初始化结果存储
results = [];

% 2. 遍历每个区域，提取并统计
for r = 1:length(regions)
    rid = regions(r);
    mask_rid = (mask == rid);  % 当前区域的掩膜（逻辑矩阵）
    
        for t = 1:5
            % 提取当前碳库、当前时间点、当前区域的所有像元值
            data = kb(:,:,t);
            values = data(mask_rid);
            values = values(~isnan(values)); % 去除 NaN
            q1 = prctile(values,25);
            q3 = prctile(values,75);
            iqr = q3 - q1;
            lower = q1 - 1.5*iqr;
            upper = q3 + 1.5*iqr;
            idx_ok = (values>=lower)&(values<=upper);
            values=values(idx_ok);            
            if isempty(values)
                mean_val = NaN; std_val = NaN; p05_val = NaN; p95_val = NaN;
            else
                mean_val = mean(values);
                std_val = std(values);
                p05_val = prctile(values,5);
                p95_val = prctile(values,95);
            end
            
            results = [results;
                {region_names{r}, carbon_vars{1}, time_env(t), ... % time_labels{t}, ...
                mean_val, std_val, p05_val, p95_val}];
        end
        
        for t = 1:5
            % 提取当前碳库、当前时间点、当前区域的所有像元值
            data = kl(:,:,t);
            values = data(mask_rid);
            values = values(~isnan(values)); % 去除 NaN
            q1 = prctile(values,25);
            q3 = prctile(values,75);
            iqr = q3 - q1;
            lower = q1 - 1.5*iqr;
            upper = q3 + 1.5*iqr;
            idx_ok = (values>=lower)&(values<=upper);
            values=values(idx_ok);            
            if isempty(values)
                mean_val = NaN; std_val = NaN; p05_val = NaN; p95_val = NaN;
            else
                mean_val = mean(values);
                std_val = std(values);
                p05_val = prctile(values,5);
                p95_val = prctile(values,95);
            end
            
            results = [results;
                {region_names{r}, carbon_vars{2}, time_env(t), ... % time_labels{t}, ...
                mean_val, std_val, p05_val, p95_val}];
        end
        
        for t = 1:5
            % 提取当前碳库、当前时间点、当前区域的所有像元值
            data = kp(:,:,t);
            values = data(mask_rid);
            values = values(~isnan(values)); % 去除 NaN
            q1 = prctile(values,25);
            q3 = prctile(values,75);
            iqr = q3 - q1;
            lower = q1 - 1.5*iqr;
            upper = q3 + 1.5*iqr;
            idx_ok = (values>=lower)&(values<=upper);
            values=values(idx_ok);            
            if isempty(values)
                mean_val = NaN; std_val = NaN; p05_val = NaN; p95_val = NaN;
            else
                mean_val = mean(values);
                std_val = std(values);
                p05_val = prctile(values,5);
                p95_val = prctile(values,95);
            end
            
            results = [results;
                {region_names{r}, carbon_vars{3}, time_env(t), ... % time_labels{t}, ...
                mean_val, std_val, p05_val, p95_val}];
        end
        
        for t = 1:5
            % 提取当前碳库、当前时间点、当前区域的所有像元值
            data = kq(:,:,t);
            values = data(mask_rid);
            values = values(~isnan(values)); % 去除 NaN
            
            q1 = prctile(values,25);
            q3 = prctile(values,75);
            iqr = q3 - q1;
            lower = q1 - 1.5*iqr;
            upper = q3 + 1.5*iqr;
            idx_ok = (values>=lower)&(values<=upper);
            values=values(idx_ok);
            if isempty(values)
                mean_val = NaN; std_val = NaN; p05_val = NaN; p95_val = NaN;
            else
                mean_val = mean(values);
                std_val = std(values);
                p05_val = prctile(values,5);
                p95_val = prctile(values,95);
            end
            
            results = [results;
                {region_names{r}, carbon_vars{4}, time_env(t), ... % time_labels{t}, ...
                mean_val, std_val, p05_val, p95_val}];
        end
        
        for t = 1:5
            % 提取当前碳库、当前时间点、当前区域的所有像元值
            data = cue(:,:,t);
            values = data(mask_rid);
            values = values(~isnan(values)); % 去除 NaN
            q1 = prctile(values,25);
            q3 = prctile(values,75);
            iqr = q3 - q1;
            lower = q1 - 1.5*iqr;
            upper = q3 + 1.5*iqr;
            idx_ok = (values>=lower)&(values<=upper);
            values=values(idx_ok);
            if isempty(values)
                mean_val = NaN; std_val = NaN; p05_val = NaN; p95_val = NaN;
            else
                mean_val = mean(values);
                std_val = std(values);
                p05_val = prctile(values,5);
                p95_val = prctile(values,95);
            end
            
            results = [results;
                {region_names{r}, carbon_vars{5}, time_env(t), ... % time_labels{t}, ...
                mean_val, std_val, p05_val, p95_val}];
        end
end

save('H:\数据库\碳源汇预测曲线3\典型区\分区统计参数.mat', 'results', '-v6');

% ============================================
% 计算不同增温时间的各过程通量
% ============================================
% 前面的5个参数值保留继续用
kb = single(zeros([1837283,5]));
kl = single(zeros([1837283,5]));
kp = single(zeros([1837283,5]));
kq = single(zeros([1837283,5]));
cue = single(zeros([1837283,5]));
kb(:,1) = param_hist(:,1);
kl(:,1) = param_hist(:,2);
kp(:,1) = param_hist(:,3);
kq(:,1) = param_hist(:,4);
cue(:,1) = param_hist(:,5);
kb(:,2:5) = param_SSP(:,1,:);
kl(:,2:5) = param_SSP(:,2,:);
kp(:,2:5) = param_SSP(:,3,:);
kq(:,2:5) = param_SSP(:,4,:);
cue(:,2:5) = param_SSP(:,5,:);
clearvars -except kb kl kp kq cue

% 只有不是变NPP就用这个碳输入
load('H:\数据库\随机森林相对重要性\全球所有点的历史环境值.mat')
NPP_HIST = new_ENVI_hist(:,7);
lat= new_ENVI_hist(:,1);
lon= new_ENVI_hist(:,2);
clear new_ENVI_hist

% Qmax %  g C kg-1  用已知的栅格图，需要转为和参数一样的数据形式，一维
% Qmax 这个在原始NC数据（H:\BaiduSyncdisk\基金\青年基金\来自文献的残体汇编数据\MOC数据\Global stocks and capacity of mineral-associated soil organic carbon\global_MOC.nc）
% 0.5度 由原来的数据直接进行单位转换，别的每没变
[Qmax0_30, R] = geotiffread('H:\数据库\Qmax0-30cm.tif');
[row, col] = geographicToDiscrete(R, lat,lon);
Qmax0_30 = Qmax0_30(sub2ind(size(Qmax0_30), row, col)); %取出对应位置的值

% 初始碳库  用之前输入的各碳库的值，需要转为和参数一样的数据形式，一维
% DOC  这个没变 是使用重采样后的BD 单位转换后的结果
[DOC, R] = geotiffread('H:\数据库\估计的DOC0-30cm_30-100cm.tif');   % 0-30cm 30-100cm 两层  
[row, col] = geographicToDiscrete(R, lat,lon); % 这个超出范围了
row(isnan(row)) = 1;
col(isnan(col)) = 1;
temp=DOC(:,:,1);
DOC30 = temp(sub2ind(size(temp), row, col)); %取出对应位置的值
clear DOC

% MAOC  这个没变  是使用重采样后的BD 单位转换后的结果
[MAOC0_30, R] = geotiffread('H:\数据库\MAOC0-30cm.tif');  
[row, col] = geographicToDiscrete(R, lat,lon);
MAOC0_30 = MAOC0_30(sub2ind(size(MAOC0_30), row, col));

% POC 这个变了 得将20cm间隔的做平均 
[POC, R] = geotiffread('H:\数据库\估计的POC0-100cm每隔20cm一层.tif');   % 0-100cm 以20cm间隔
POC0_30 = mean(cat(3,POC(:,:,1),POC(:,:,1),POC(:,:,2)), 3,'omitnan');
clear POC
[row, col] = geographicToDiscrete(R, lat,lon); % 这个超出范围了
row(isnan(row)) = 1;
col(isnan(col)) = 1;
POC30 = POC0_30(sub2ind(size(POC0_30), row, col)); %取出对应位置的值
clear POC0_30

% MBC 这个没变
% 只按0-30  30-100cm分层 单位 mgC kg-1Soil;
[mbc0_30, R] = geotiffread( 'J:\数据备份-2023.12.18\global pattern\MBC30cm.tif');
mask2 = (mbc0_30 <= 0);
mbc0_30(mask2) = NaN;
mbc0_30 = mbc0_30/1000; % gC kg-1Soil;
clear mask2
[row, col] = geographicToDiscrete(R, lat,lon);
row(isnan(row)) = 1;
col(isnan(col)) = 1;
mbc30 = mbc0_30(sub2ind(size(mbc0_30), row, col));
clear mbc0_30

% 合并初始值
wtchs30 = [POC30,DOC30,mbc30,MAOC0_30];
clear POC30 DOC30 mbc30 MAOC0_30

indir='H:\数据库\随机森林相对重要性\未来情景的预测值\';
yt = {'2140','4160','6180','8100'};
NPP_SSP = single(zeros(1837283,4));
idh = 1;
for jj=1:4
	files = [indir, 'SSP585_',yt{jj},'的环境值.mat'];
	temp = load(files);
	fieldNames = fieldnames(temp);
	originalVarName = fieldNames{1};
	temp = temp.(originalVarName);
	NPP_SSP(:,idh)=temp(:,7);
	idh = idh+1;
end

in0_30 = single(zeros(1837283,5));
in0_30(:,1) = 3.474300E-03*NPP_HIST + 4.682669E-04;   % 全变和只变NPP的时候用这个数据
in0_30(:,2:5) = 3.474300E-03*NPP_SSP + 4.682669E-04;   % 全变和只变NPP的时候用这个数据
clear NPP_HIST NPP_SSP

% 常数值
pi=0.66;
kml = 0.001; % g C kg-1 day-1
kml = repmat(kml, 1837283, 1);
    
load('H:\数据库\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量全变year.mat')

target_years = [2026, 2039, 2060, 2078, 2094];
Global_cpool_pre30 = Global_cpool_pre30(:,:,target_years-2020);
Global_cpool_pre30 = Global_cpool_pre30+wtchs30;

flux = single(zeros(1837283,7,6)); % 7个通量，6个时间点

for i =1:6
    if i == 1       
        % 历史碳库
        P = wtchs30(:,1);
        L = wtchs30(:,2);
        B = wtchs30(:,3);
        M = wtchs30(:,4);
        tt=1;
    end
    if i == 2 || i == 3 
        % 未来碳库
        P = Global_cpool_pre30(:,1,i-1);
        L = Global_cpool_pre30(:,2,i-1);
        B = Global_cpool_pre30(:,3,i-1);
        M = Global_cpool_pre30(:,4,i-1);
        tt=2;
    end
    if i == 4 
        % 未来碳库
        P = Global_cpool_pre30(:,1,i-1);
        L = Global_cpool_pre30(:,2,i-1);
        B = Global_cpool_pre30(:,3,i-1);
        M = Global_cpool_pre30(:,4,i-1);
        tt=3;
    end
    if i == 5 
        % 未来碳库
        P = Global_cpool_pre30(:,1,i-1);
        L = Global_cpool_pre30(:,2,i-1);
        B = Global_cpool_pre30(:,3,i-1);
        M = Global_cpool_pre30(:,4,i-1);
        tt=4;
    end
	if i == 6 
        % 未来碳库
        P = Global_cpool_pre30(:,1,i-1);
        L = Global_cpool_pre30(:,2,i-1);
        B = Global_cpool_pre30(:,3,i-1);
        M = Global_cpool_pre30(:,4,i-1);
        tt=5;
    end
    % 向量化计算通量（.*表示元素乘法）
    klm = kml.*kq(:,tt);  % 元素乘法
    Fpl = kp(:,tt).*P;
    Fbm = kb(:,tt).*B;   % ,注意这里的kbd为kbd*mbc
    Fml = kml.*M./Qmax0_30;
    Flm = klm.*L.*(1-M./Qmax0_30);
    Flb = kl(:,tt).*L;
    CO2 = (1-cue(:,tt)).*Flb;
    flux(:,:,i) = [in0_30(:,tt),CO2,Fpl,Flb,Fbm,Flm,Fml];
end
clearvars -except flux

I = single(zeros([2160,4320,6]));
CO2 = single(zeros([2160,4320,6]));
Fpl = single(zeros([2160,4320,6]));
Flb = single(zeros([2160,4320,6]));
Fbm = single(zeros([2160,4320,6]));
Flm = single(zeros([2160,4320,6]));
Fml = single(zeros([2160,4320,6]));
I(:) = NaN;
CO2(:) = NaN;
Fpl(:) = NaN;
Flb(:) = NaN;
Fbm(:) = NaN;
Flm(:) = NaN;
Fml(:) = NaN;

indir='H:\数据库\随机森林相对重要性\未来情景的预测值\';
load([indir, 'SSP126_2140的环境值.mat']);
for ii=1:length(envi_SSP126_2140)
    I(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,:) = flux(ii,1,:);
    CO2(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,:) = flux(ii,2,:);
    Fpl(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,:) = flux(ii,3,:);
    Flb(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,:) = flux(ii,4,:);
    Fbm(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,:) = flux(ii,5,:);
    Flm(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,:) = flux(ii,6,:);
    Fml(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,:) = flux(ii,7,:);
end

clearvars -except I CO2 Fpl Flb Fbm Flm Fml
sz = size(CO2);
FLUX = zeros(sz(1), sz(2), sz(3), 7, 'single');
FLUX(:,:,:,1) = I;
FLUX(:,:,:,2) = CO2;
FLUX(:,:,:,3) = Fpl;
FLUX(:,:,:,4) = Flb;
FLUX(:,:,:,5) = Fbm;
FLUX(:,:,:,6) = Flm;
FLUX(:,:,:,7) = Fml;
clearvars -except FLUX

% 典型区掩膜======================
[mask, ~] = geotiffread('H:\数据库\碳源汇预测曲线3\典型区\MASK2.tif');
 % 1 亚马孙  2 北方森林  3 北方森林北美  4 刚果  5 北极苔原 6 北极苔原北美
% 合并相同气候区
mask(mask==4)=1; % 1 亚马孙 4 刚果  ==热带雨林
mask(mask==3)=2; % 2 北方森林 3 北方森林北美  ==北方森林
mask(mask==5)=3; % 5 北极苔原 6 北极苔原北美  ==北极苔原
mask(mask==6)=3; % 5 北极苔原 6 北极苔原北美  ==北极苔原
 
% 1. 定义参数和区域列表
% 区域编号（按你说的顺序）
% regions = [1 2 3 4 5 6];
regions = [1 2 3]; % 合并相同气候区
% region_names = {'Amazon','Boreal_Eurasia','Boreal_NorthAmerica','Congo','Tundra_Eurasia','Tundra_NorthAmerica'};
region_names = {'Tropical_rainforest','Boreal_Eurasia','Tundra_Eurasia'}; % 合并相同气候区

% 变量名
carbon_vars = {'I', 'CO2', 'Fpl', 'Flb', 'Fbm', 'Flm', 'Fml'};

% 时间点标签（对应第4维/第3维）
time_env = {'hist','2026', '2039', '2060', '2078', '2094'};

% 初始化结果存储
results = [];

% 2. 遍历每个区域，提取并统计
for r = 1:length(regions)
    rid = regions(r);
    mask_rid = (mask == rid);  % 当前区域的掩膜（逻辑矩阵）
    for t = 1:length(time_env)
        for v = 1:length(carbon_vars)
            % 提取当前碳库、当前时间点、当前区域的所有像元值
            data = FLUX(:,:,t,v);
            values = data(mask_rid);
            values = values(~isnan(values)); % 去除 NaN
            q1 = prctile(values,25);
            q3 = prctile(values,75);
            iqr = q3 - q1;
            lower = q1 - 1.5*iqr;
            upper = q3 + 1.5*iqr;
            idx_ok = (values>=lower)&(values<=upper);
            values=values(idx_ok);            
            if isempty(values)
                mean_val = NaN; std_val = NaN; p05_val = NaN; p95_val = NaN;
            else
                mean_val = mean(values);
                std_val = std(values);
                p05_val = prctile(values,5);
                p95_val = prctile(values,95);
            end
            
            results = [results;
                {region_names{r}, time_env(t), carbon_vars{v},... % time_labels{t}, ...
                mean_val, std_val, p05_val, p95_val}];
        end
    end
end

save('H:\数据库\碳源汇预测曲线3\典型区\分区统计flux.mat', 'results', '-v6');
     

%% 将预测一维的计算结果转为二维图，然后在提取对应特征区域的平均值变化
% 修改下面的部分
% 历史环境值===========================
indir_env = 'F:\';
indir_ssp = 'H:\数据库\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\';
indir_npp = 'H:\数据库\NPP_NDVI_MODIS\Global_Mean_NPP_2003_2023_MODIS重采样0.0833.tif';

% 1.1 加载ET0/AI/MAT/MAP
load([indir_env, '参考蒸散量ET0数据\gaez-v5-data\历史数据\1980-2022年球ET0平均值.mat']); % meanET0
load([indir_env, '干旱指数\用MAP计算的AI\1980-2022年球AI平均值.mat']); % meanAI1
load([indir_env, 'ERA5气候数据\20251111下载的新的\1980-2024年球年平均气温平均值.mat']); % meanMAT
load([indir_env, 'ERA5气候数据\20251111下载的新的\1980-2024年球年降水量平均值.mat']); % meanMAP

% 1.2 读取NPP
[meanNPP, ~] = geotiffread(indir_npp);

% 1.3 统一预处理MAT/MAP（维度匹配）
% MAT预处理
result=meanMAT;
result = result(:,1:1800,:);
result = rot90(result, -1); 
result = flip(result, 2);
result= [result(:,1801:3600,:),result(:,1:1800,:)]; 
result2 = imresize(result, [2160 4320], 'nearest');
meanMAT=result2; clear result result2;

% MAP预处理
result=meanMAP;
result = result(:,1:1800,:);
result = rot90(result, -1); 
result = flip(result, 2);
result= [result(:,1801:3600,:),result(:,1:1800,:)]; 
result2 = imresize(result, [2160 4320], 'nearest');
meanMAP=result2; clear result result2;

mask=(~isnan(meanMAT)) & (meanMAT > -14) & (meanMAT < 29) & ...
        (~isnan(meanMAP)) & (meanMAP > 200) & (meanMAP < 3000) & (~isnan(meanAI1)) & ...
        (meanAI1 > 600) & (meanAI1 < 29000) & (meanET0 ~= -9) & (meanET0 > 160) & (meanET0 < 2800) & ...
        (meanNPP ~= 0) & (~isnan(meanNPP)) & (meanNPP > 0.006) & (meanNPP < 1.57);

meanET0=double(meanET0);
meanAI1(~mask)=NaN;
meanET0(~mask)=NaN;
meanMAP(~mask)=NaN;
meanMAT(~mask)=NaN;
meanNPP(~mask)=NaN;


% 未来情景的环境值========================
indir='H:\数据库\随机森林相对重要性\未来情景的预测值\';
NPP_SSP = single(zeros(1837283,4));
AI_SSP = single(zeros(1837283,4));
ET_SSP = single(zeros(1837283,4));
MAP_SSP = single(zeros(1837283,4));
MAT_SSP = single(zeros(1837283,4));
yt = {'2140','4160','6180','8100'};
idh = 1;
for jj=1:4
	files = [indir, 'SSP585_',yt{jj},'的环境值.mat'];
	temp = load(files);
	fieldNames = fieldnames(temp);
	originalVarName = fieldNames{1};
	temp = temp.(originalVarName);
	NPP_SSP(:,idh)=temp(:,7);
    AI_SSP(:,idh)=temp(:,6);
    ET_SSP(:,idh)=temp(:,5);
    MAP_SSP(:,idh)=temp(:,4);
    MAT_SSP(:,idh)=temp(:,3);
	idh = idh+1;
end

% 未来的预测碳库值====================================
indir='H:\数据库\';
load([indir, 'SSP126_2140的环境值.mat']);
info = geotiffinfo('F:\参考蒸散量ET0数据\gaez-v5-data\temp.tif');
[temp, R] = geotiffread('F:\参考蒸散量ET0数据\gaez-v5-data\temp.tif');

% 20260524
indir2='H:\数据库\全球只变单一变量预测的参数4\';
load([indir2, '全球碳库SSP预测0-30增量全变year.mat']); % 点 4个库 时间序列 

indir='H:\数据库\增量碳库预测栅格4\'; 
yt = {'1.5C','2.0C','3.0C','4.0C','5.0C','2100y'};

% 碳库数据与环境数据
carbon30 = single(zeros([2160,4320,5,80]));

for kk=1:80 %kk=1:6
    param_SSP_tif30 = single(zeros([size(temp),4])); % 4个库
    param_SSP_tif30(:) = NaN;
    for ii=1:length(envi_SSP126_2140)
        param_SSP_tif30(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,:) = Global_cpool_pre30(ii,:,kk); % 2140 ssp126
    end
    param_SSP_tif30SOC = sum(param_SSP_tif30,3);
    
    carbon30(:,:,1:4,kk) = param_SSP_tif30;
    carbon30(:,:,5,kk) = param_SSP_tif30SOC;
end
clear param_SSP_tif30 param_SSP_tif30SOC Global_cpool_pre30
clear temp

MAT = single(zeros([2160,4320,5]));
MAP = single(zeros([2160,4320,5]));
ET = single(zeros([2160,4320,5]));
AI = single(zeros([2160,4320,5]));
NPP = single(zeros([2160,4320,5]));
MAT(:) = NaN;
MAP(:) = NaN;
ET(:) = NaN;
AI(:) = NaN;
NPP(:) = NaN;
for ii=1:length(envi_SSP126_2140)
    MAT(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,2:5) = MAT_SSP(ii,:);
    MAP(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,2:5) = MAP_SSP(ii,:);
    ET(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,2:5) = ET_SSP(ii,:);
    AI(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,2:5) = AI_SSP(ii,:);
    NPP(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,2:5) = NPP_SSP(ii,:);
end

MAT(:,:,1) = meanMAT;
MAP(:,:,1) = meanMAP;
ET(:,:,1) = meanET0;
AI(:,:,1) = meanAI1;
NPP(:,:,1) = meanNPP;

% 只保留 ss、dd，其他变量清空
clearvars -except MAT MAP ET AI NPP carbon30 carbon100


% 典型区掩膜======================
[mask, ~] = geotiffread('H:\数据库\碳源汇预测曲线3\典型区\MASK2.tif');
 % 1 亚马孙  2 北方森林  3 北方森林北美  4 刚果  5 北极苔原 6 北极苔原北美
% 合并相同气候区
mask(mask==4)=1; % 1 亚马孙 4 刚果  ==热带雨林
mask(mask==3)=2; % 2 北方森林 3 北方森林北美  ==北方森林
mask(mask==5)=3; % 5 北极苔原 6 北极苔原北美  ==北极苔原
mask(mask==6)=3; % 5 北极苔原 6 北极苔原北美  ==北极苔原


% 新加的一个掩膜 去除苔原和一些土地类型 2026 09 16
[TOC301_5, ~] = geotiffread('H:\数据库\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_1.5C增量TOC.tif');
[EAS, ~] = geotiffread('H:\数据库\ESA_CCI_LandCoverData\重采样.tif');  % 掩膜掉10 11 12 20 30 190 200 201 202 210 220
[KOP, ~] = geotiffread('H:\数据库\气候类型数据\koppen_geiger_tif\1991_2020\重采样.tif');  % 掩膜掉29
% 条件1：TOC301_5为NaN
cond1 = isnan(TOC301_5);
% 条件2：EAS属于指定类别
eas_list = [10,11,12,20,30,190,200,201,202,210,220];
cond2 = ismember(EAS, eas_list);
% 条件3：KOP等于29
cond3 = (KOP == 29);
% 三者并集：任意条件成立则掩膜位置为0
mask_union = cond1 | cond2 | cond3;
% 生成最终掩膜：mask=0为剔除，mask=1为保留
mask3 = mask < 4;
mask3 = ~mask_union.*mask3;
mask = mask .* uint8(mask3);
mask(mask==0)=15;


% 1. 定义参数和区域列表
% 区域编号（按你说的顺序）
% regions = [1 2 3 4 5 6];
regions = [1 2 3]; % 合并相同气候区
% region_names = {'Amazon','Boreal_Eurasia','Boreal_NorthAmerica','Congo','Tundra_Eurasia','Tundra_NorthAmerica'};
region_names = {'Tropical_rainforest','Boreal_Eurasia','Tundra_Eurasia'}; % 合并相同气候区

% 碳库变量名（对应第3维：POC DOC MBC MAOC TOC）
carbon_vars = {'POC','DOC','MBC','MAOC','TOC'};

% 环境变量
% env_vars = {'MAT','MAP','ET','AI','NPP'};
env_vars = {'MAT','NPP'}; %20260916

% 时间点标签（对应第4维/第3维）
% time_labels = {'1.5C','2.0C','3.0C','4.0C','5.0C'};
time_labels = (2021:2100);
time_env = {'hist','2020-2040','2040-2060','2060-2080','2080-2100'};
time_env2 = [2010,2040,2060,2080,2100];
% 20260916
% time_labels = (2061:2100);
% time_env = {'hist','2020-2040','2040-2060'};
% time_env2 = [2010,2040,2060];

% 初始化结果存储
results = [];

% 2. 遍历每个区域，提取并统计
for r = 1:length(regions)
    rid = regions(r);
    mask_rid = (mask == rid);  % 当前区域的掩膜（逻辑矩阵）
    
    % --------------------------
    % 处理碳库 carbon30（0-30cm）
    % --------------------------
    for c = 1:length(carbon_vars)
        results = [results;
            {region_names{r}, 'carbon30', carbon_vars{c}, 2010, ... % time_labels{t}, ...
            0, 0, 0, 0}];
        for t = 1:length(time_labels)
            % 提取当前碳库、当前时间点、当前区域的所有像元值
            data = carbon30(:,:,c,t);
            values = data(mask_rid);
            values = values(~isnan(values)); % 去除 NaN
            
            if isempty(values)
                mean_val = NaN; std_val = NaN; p05_val = NaN; p95_val = NaN;
            else
                % 先去异常值
                mu0 = mean(values);
                sigma0 = std(values);
                lower_thr = mu0 - 3*sigma0;
                upper_thr = mu0 + 3*sigma0;
                vec_clean = values(values >= lower_thr & values <= upper_thr);
                n_after = length(vec_clean);
                % 求平均和95置信区间
                mean_val = mean(vec_clean);
                std_val = std(vec_clean);
                se = std_val / sqrt(n_after);
                t_crit = tinv(0.975, n_after-1);
                p05_val = mean_val - t_crit*se;
                p95_val = mean_val + t_crit*se;
            end
            
            results = [results;
                {region_names{r}, 'carbon30', carbon_vars{c}, time_labels(t), ... % time_labels{t}, ...
                mean_val, std_val, p05_val, p95_val}];
        end
    end
    
    
    % --------------------------
    % 处理环境变量 MAT/MAP/ET/AI/NPP
    % --------------------------
% %     env_data_list = {MAT, MAP, ET, AI, NPP};
%     env_data_list = {MAT, NPP}; % 20260916  环境和碳库分别计算
%     for e = 1:length(env_vars)
%         env_data = env_data_list{e}; % [2160,4320,5]
%         for t = 1:size(env_data,3)
%             data = env_data(:,:,t);
%             values = data(mask_rid);
%             values = values(~isnan(values));
%             
%             if isempty(values)
%                 mean_val = NaN; std_val = NaN; p05_val = NaN; p95_val = NaN;
%             else
%                 % 先去异常值
%                 mu0 = mean(values);
%                 sigma0 = std(values);
%                 lower_thr = mu0 - 3*sigma0;
%                 upper_thr = mu0 + 3*sigma0;
%                 vec_clean = values(values >= lower_thr & values <= upper_thr);
%                 n_after = length(vec_clean);
%                 % 求平均和95置信区间
%                 mean_val = mean(vec_clean);
%                 std_val = std(vec_clean);
%                 se = std_val / sqrt(n_after);
%                 t_crit = tinv(0.975, n_after-1);
%                 p05_val = mean_val - t_crit*se;
%                 p95_val = mean_val + t_crit*se;
%             end
%             
%             results = [results;
%                 {region_names{r}, env_vars{e}, time_env{t}, time_env2(t), ...  %time_labels{t}, ...
%                 mean_val, std_val, p05_val, p95_val}];
%         end
%     end
end

save('H:\数据库\碳源汇预测曲线3\典型区\分区统计year_hebing4.mat', 'results', '-v6');

% 3. 把结果转为表格并保存
T = cell2table(results, ...
    'VariableNames', {'Region','VariableType','CarbonPool','TimePoint', ...
    'Mean','Std','P05','P95'});

% 保存到 CSV，方便后续在 R 中绘图
writetable(T, 'H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\碳源汇预测曲线3\典型区\分区统计year_hebing4.csv');

% 将数据整理为origin 输入的形式
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\碳源汇预测曲线3\典型区\分区统计year_hebing.mat')
T = cell2table(results, ...
    'VariableNames', {'Region','VariableType','CarbonPool','Year', ...
    'Mean','Std','P05','P95'});

% 转成数值/字符数组
RegionList     = unique(T.Region);
CarbonTypeList = {'carbon30','carbon100'};
VarList        = {'POC','DOC','MBC','MAOC','TOC','MAT','MAP','ET','AI','NPP'};
nRegion        = length(RegionList);
nCarbonType    = length(CarbonTypeList);
nVars          = length(VarList);
allYears       = unique(T.Year);
nYear          = length(allYears);

% 环境变量的目标年份（对应你提供的环境数据）
envYears = [2010, 2040, 2060, 2080, 2100];

% ====================== 2. 按区域+碳库深度分组，整理均值和标准差 ======================
for iR = 1:nRegion
    regionName = RegionList{iR};
    for iC = 1:nCarbonType
        carbonType = CarbonTypeList{iC};
        
        % 筛选当前区域 + 当前碳库深度的数据
        idx = strcmp(T.Region, regionName) & ...
              strcmp(T.VariableType, carbonType);
        T_sub = T(idx,:);
        
        % 初始化当前组的矩阵：Year x 变量（第一列是年份）
        data_mean = [allYears, nan(nYear, nVars)];
        data_std  = [allYears, nan(nYear, nVars)];
        
        % --------------------------
        % 填充碳库数据（POC/LMWC/MBC/MAOC/TOC）
        % --------------------------
        for iVar = 1:5
            varName = VarList{iVar};
            idx_var = strcmp(T_sub.CarbonPool, varName);
            T_var = T_sub(idx_var,:);
            
            if ~isempty(T_var)
                [~, loc] = ismember(T_var.Year, allYears);
                data_mean(loc, iVar+1) = T_var.Mean;
                data_std(loc, iVar+1)  = T_var.Std;
            end
        end
        
        % --------------------------
        % 填充环境数据（MAT/MAP/ET/AI/NPP）
        % --------------------------
        for iVar = 6:nVars
            varName = VarList{iVar};
            idx_var = strcmp(T.VariableType, varName) & strcmp(T.Region, regionName);
            T_var = T(idx_var,:);
            
            if ~isempty(T_var)
                % 只保留环境变量对应的年份
                idx_env = ismember(T_var.Year, envYears);
                T_var_env = T_var(idx_env,:);
                
                [~, loc] = ismember(T_var_env.Year, allYears);
                data_mean(loc, iVar+1) = T_var_env.Mean;
                data_std(loc, iVar+1)  = T_var_env.Std;
            end
        end
        
        % 构造变量名，比如 Amazon_carbon30_mean, Amazon_carbon30_std
        baseName = [regionName '_' carbonType];
        eval([baseName '_mean = data_mean;']);
        eval([baseName '_std  = data_std;']);
    end
end
 
 
        
%% 将预测一维的计算结果转为二维图
% 修改BaoHe3 的 4部分
indir='H:\数据库\随机森林相对重要性\未来情景的预测值\';
load([indir, 'SSP126_2140的环境值.mat']);
info = geotiffinfo('F:\参考蒸散量ET0数据\gaez-v5-data\temp.tif');
[temp, R] = geotiffread('F:\参考蒸散量ET0数据\gaez-v5-data\temp.tif');

% 20260425 修改
indir2='H:\数据库\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\';
load([indir2, '全球碳库SSP预测0-30增量全变.mat']); % 点 4个库 6个时间点 

indir='H:\数据库\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\'; 
yt = {'1.5C','2.0C','3.0C','4.0C','5.0C','2100y'};

    for kk=1:6
        param_SSP_tif30 = single(zeros([size(temp),4])); % 4个库 
        param_SSP_tif30(:) = NaN;
        for ii=1:length(envi_SSP126_2140)
            param_SSP_tif30(envi_SSP126_2140(ii,1)+1,envi_SSP126_2140(ii,2)+1,:) = Global_cpool_pre30(ii,:,kk); % 2140 ssp126
        end
        geotiffwrite([indir, '全球碳库SSP预测30_',yt{kk},'增量.tif'],param_SSP_tif30,R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag)
        param_SSP_tif30SOC = sum(param_SSP_tif30,3);
        geotiffwrite([indir, '全球碳库SSP预测30_',yt{kk},'增量TOC.tif'],param_SSP_tif30SOC,R, 'GeoKeyDirectoryTag', info.GeoTIFFTags.GeoKeyDirectoryTag)
    end

% 取出随机点的值===========
indir='H:\数据库\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格3\';
load('H:\数据库\随机森林相对重要性\新选的更多点的全球参数和环境值_已去除异常值.mat')
latlon = cleaned_new_param_envi(:,1:2);
temp3 = single(zeros([285199,5]));
% [temp, R] = geotiffread([indir, '全球碳库SSP预测30_SSP585_8100增量.tif']);
% [temp, R] = geotiffread([indir, '全球碳库SSP预测30_SSP370_8100增量.tif']);
[temp, R] = geotiffread([indir, '全球碳库SSP预测30_SSP126_8100增量.tif']);
[row, col] = geographicToDiscrete(R, latlon(:,1),latlon(:,2));
temp2 = temp(:,:,1);
temp3(:,1) = temp2(sub2ind(size(temp2), row, col)); %取出对应位置的值
temp2 = temp(:,:,2);
temp3(:,2) = temp2(sub2ind(size(temp2), row, col)); %取出对应位置的值
temp2 = temp(:,:,3);
temp3(:,3) = temp2(sub2ind(size(temp2), row, col)); %取出对应位置的值
temp2 = temp(:,:,4);
temp3(:,4) = temp2(sub2ind(size(temp2), row, col)); %取出对应位置的值

% [temp, ~] = geotiffread([indir, '全球碳库SSP预测30_SSP585_8100增量TOC.tif']);
% [temp, ~] = geotiffread([indir, '全球碳库SSP预测30_SSP370_8100增量TOC.tif']);
[temp, ~] = geotiffread([indir, '全球碳库SSP预测30_SSP126_8100增量TOC.tif']);
temp3(:,5) = temp(sub2ind(size(temp), row, col)); %取出对应位置的值

temp3 = [temp3,cleaned_new_param_envi(:,16:20)];
nan_mask = any(isnan(temp3), 2);
temp3=temp3(~nan_mask,:);
% save([indir, '新选的更多点的_全球碳库SSP预测30_SSP585_8100增量.mat'], 'temp3', '-v6');
% save([indir, '新选的更多点的_全球碳库SSP预测30_SSP370_8100增量.mat'], 'temp3', '-v6');
save([indir, '新选的更多点的_全球碳库SSP预测30_SSP126_8100增量.mat'], 'temp3', '-v6');



%% 2.计算SSP585在不同时间点的结果,这个与下面的差别在于是否用数据集的值作为初始值。
% 下面的那个是使用的稳态值作为初始值
% 修改BaoHe3 的 3.2 部分，不区分SSP
% 3.2、将预测的未来的全球SSP的参数值代入模型 预测未来碳库  以及与历史的差值
% 修改第3部分，因为第2部分用的随机选的点代表性较差，这里计算的所有像元的
% 每个点一个曲线，不存储曲线值，只存储每个时间段的最后一张图，碳库数据 ==============================
% 读取所有数据
indir='H:\数据库\随机森林相对重要性\未来情景的预测值\';
indir2='H:\数据库\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数\';  % 只变一个变量预测值
param_name = {'kbd_0_30_83', 'kbd_30_100_83', 'kl_0_30_83', 'kl_30_100_83', 'kp_0_30_83', 'kp_30_100_83', 'kq_0_30_83', 'kq_30_100_83', 'CUE'};
yt = {'2140','4160','6180','8100'};
hh = 1837283;  %全变行数
% hh = 1581186; % 只变一个变量预测值  行数 MAT MAP ET
% hh = 1582872; % 只变一个变量预测值  行数 AI 
% hh = 1833546; % 只变一个变量预测值  行数   NPP
param_SSP = single(zeros(hh,9,4)); % 只变一个变量预测值  9个变量 4个时间点
for ii =1:9
    for jj=1:4
        files = [indir, 'Global_',param_name{ii},'_SSP585_',yt{jj},'.mat'];
%         files = [indir2, 'Global_',param_name{ii},'_SSP585_',yt{jj},'qu_MAT.mat']; % 只变一个变量预测值
%         files = [indir2, 'Global_',param_name{ii},'_SSP585_',yt{jj},'qu_MAP.mat']; % 只变一个变量预测值
%         files = [indir2, 'Global_',param_name{ii},'_SSP585_',yt{jj},'qu_ET.mat']; % 只变一个变量预测值
%         files = [indir2, 'Global_',param_name{ii},'_SSP585_',yt{jj},'qu_AI.mat']; % 只变一个变量预测值
%         files = [indir2, 'Global_',param_name{ii},'_SSP585_',yt{jj},'qu_NPP.mat']; % 只变一个变量预测值
        temp = load(files);
        fieldNames = fieldnames(temp);
        originalVarName = fieldNames{1};
        temp2 = temp.(originalVarName);
        param_SSP(:,ii,jj)=temp2(:,3);
    end
end

NPP_SSP = single(zeros(1837283,4));
idh = 1;
for jj=1:4
	files = [indir, 'SSP585_',yt{jj},'的环境值.mat'];
	temp = load(files);
	fieldNames = fieldnames(temp);
	originalVarName = fieldNames{1};
	temp = temp.(originalVarName);
	NPP_SSP(:,idh)=temp(:,7);
	idh = idh+1;
end


% 只变一个变量预测值 取出对应的行=============
[c, ia, ib] = intersect(temp(:, 1:2), temp2(:, 1:2), 'rows');
NPP_SSP=NPP_SSP(ia,:);
temp=temp2;
% 只有不是变NPP就用这个碳输入
load('H:\数据库\随机森林相对重要性\全球所有点的历史环境值.mat')
NPP_HIST = new_ENVI_hist(ia,7);
% 只变一个变量预测值 取出对应的行=============

% 准备其他碳库数据====
% 取出所有像元的经纬度
info = geotiffinfo('F:\参考蒸散量ET0数据\gaez-v5-data\temp.tif');
[lat,lon]=pix2latlon(info.RefMatrix,temp(:,1)+1,temp(:,2)+1); % 之前的索引从0开始的,在MATLAB中需要加1

% Qmax %  g C kg-1  用已知的栅格图，需要转为和参数一样的数据形式，一维
% 0.5度 由原来的数据直接进行单位转换，别的每没变
[Qmax0_30, R] = geotiffread('H:\数据库\Qmax0-30cm.tif');
[Qmax30_100, ~] = geotiffread('H:\数据库\Qmax30-100cm.tif');
[row, col] = geographicToDiscrete(R, lat,lon);
Qmax0_30 = Qmax0_30(sub2ind(size(Qmax0_30), row, col)); %取出对应位置的值
Qmax30_100 = Qmax30_100(sub2ind(size(Qmax30_100), row, col)); %取出对应位置的值

% 初始碳库  用之前输入的各碳库的值，需要转为和参数一样的数据形式，一维
% DOC  这个没变 是使用重采样后的BD 单位转换后的结果
[DOC, R] = geotiffread('H:\数据库\估计的DOC0-30cm_30-100cm.tif');   % 0-30cm 30-100cm 两层  
[row, col] = geographicToDiscrete(R, lat,lon); % 这个超出范围了
mask = ~(isnan(row) |  isnan(col));
temp=DOC(:,:,1);
DOC30 = single(zeros(hh,1));
DOC30(mask) = temp(sub2ind(size(temp), row(mask), col(mask))); %取出对应位置的值
DOC30(~mask) = nan;
temp=DOC(:,:,2);
DOC100 = single(zeros(hh,1));
DOC100(mask) = temp(sub2ind(size(temp), row(mask), col(mask))); 
DOC100(~mask) = nan;
clear DOC

% MAOC  这个没变  是使用重采样后的BD 单位转换后的结果
[MAOC0_30, R] = geotiffread('H:\数据库\MAOC0-30cm.tif');  
[MAOC30_100, ~] = geotiffread('H:\数据库\MAOC30-100cm.tif');
[row, col] = geographicToDiscrete(R, lat,lon);
MAOC0_30 = MAOC0_30(sub2ind(size(MAOC0_30), row, col));
MAOC30_100 = MAOC30_100(sub2ind(size(MAOC30_100), row, col));

% POC 这个变了 得将20cm间隔的做平均 
[POC, R] = geotiffread('H:\数据库\估计的POC0-100cm每隔20cm一层.tif');   % 0-100cm 以20cm间隔
POC0_30 = mean(cat(3,POC(:,:,1),POC(:,:,1),POC(:,:,2)), 3,'omitnan');
POC30_100 = mean(cat(3,POC(:,:,2),POC(:,:,3),POC(:,:,4),POC(:,:,5),POC(:,:,3),POC(:,:,4),POC(:,:,5)), 3,'omitnan');
clear POC
[row, col] = geographicToDiscrete(R, lat,lon); % 这个超出范围了
mask = ~(isnan(row) |  isnan(col));
POC30 = single(zeros(hh,1));
POC30(mask) = POC0_30(sub2ind(size(POC0_30), row(mask), col(mask))); %取出对应位置的值
POC30(~mask) = nan;
POC100 = single(zeros(hh,1));
POC100(mask) = POC30_100(sub2ind(size(POC30_100), row(mask), col(mask))); %取出对应位置的值
POC100(~mask) = nan;
clear POC0_30 POC30_100

% MBC 这个没变
% 只按0-30  30-100cm分层 单位 mgC kg-1Soil;
[mbc30_100, R] = geotiffread('J:\数据备份-2023.12.18\global pattern\MBC_30-100cm.tif');    
[mbc0_30, ~] = geotiffread( 'J:\数据备份-2023.12.18\global pattern\MBC30cm.tif');
mask2 = (mbc30_100 <= 0) | (mbc0_30 <= 0);
mbc30_100(mask2) = NaN;
mbc0_30(mask2) = NaN;
mbc30_100 = mbc30_100/1000; % gC kg-1Soil;
mbc0_30 = mbc0_30/1000; % gC kg-1Soil;
clear mask2
[row, col] = geographicToDiscrete(R, lat,lon);
mbc100 = mbc30_100(sub2ind(size(mbc30_100), row, col));
mbc30 = mbc0_30(sub2ind(size(mbc0_30), row, col));
clear mbc0_30 mbc30_100

% 合并初始值
wtchs30 = [POC30,DOC30,mbc30,MAOC0_30];
wtchs100 = [POC100,DOC100,mbc100,MAOC30_100];
clear POC30 DOC30 mbc30 MAOC0_30
clear POC100 DOC100 mbc100 MAOC30_100

% 常数值
pi=0.66;
kml = 0.001; % g C kg-1 day-1

% 运行时间 
years = 20;
Nt=years*365; % 年

% 20260524每年都有值
Global_cpool_pre30 = single(zeros(hh,4,80));   % 点 4个库 80年
Global_cpool_pre100 = single(zeros(hh,4,80));   % 点 4个库 80年

% 定义你需要提取的目标年份
% target_years = [2026, 2039, 2060, 2078, 2094, 2100];
% 20260524每年都有值
target_years = (2021:2100);

idd=1; % 截取时间点计数器
    % 每个时间段一个循环 21-40
    for tt = 1:4
        fprintf('year=%d\n',tt);
        kp30 = param_SSP(:,5,tt);
        kp100 = param_SSP(:,6,tt);
        kbd30 = param_SSP(:,1,tt);
        kbd100 = param_SSP(:,2,tt);
        kq30 = param_SSP(:,7,tt);
        kq100 = param_SSP(:,8,tt);
        kl30 = param_SSP(:,3,tt);
        kl100 = param_SSP(:,4,tt);
        cue = param_SSP(:,9,tt);
        %==========================
        in0_30 = 3.474300E-03*NPP_SSP(:,tt) + 4.682669E-04;   % 全变和只变NPP的时候用这个数据
        in30_100 = 1.342353E-03*NPP_SSP(:,tt) + 1.650560E-07;   % 全变和只变NPP的时候用这个数据
        %==========================
        % 因为只有碳输入是用NPP计算的。所以不变NPP时,碳输入也不变。
%         in0_30 = 3.474300E-03*NPP_HIST + 4.682669E-04;  % 只变一个变量时用=====================
%         in30_100 = 1.342353E-03*NPP_HIST + 1.650560E-07; % 只变一个变量时用=====================
        
        % 修改后的参数
        k_temp30=[kp30,kbd30,kq30,kl30,cue,in0_30,repmat(pi, hh, 1),repmat(kml, hh, 1),Qmax0_30,wtchs30];
        k_temp100=[kp100,kbd100,kq100,kl100,cue,in30_100,repmat(pi, hh, 1),repmat(kml, hh, 1),Qmax30_100,wtchs100];  % 30-100
        nan_mask = any(isnan(k_temp30), 2);       % 标记含NaN的行
        neg_mask = any(k_temp30< 0, 2);          % 标记含负值的行
        target_rows = nan_mask | neg_mask; % 同时满足条件的行逻辑索引
        k_temp30(target_rows, :) = NaN;
        nan_mask = any(isnan(k_temp100), 2);       % 标记含NaN的行
        neg_mask = any(k_temp100< 0, 2);          % 标记含负值的行
        target_rows = nan_mask | neg_mask; % 同时满足条件的行逻辑索引
        k_temp100(target_rows, :) = NaN;
        
        k_temp30=k_temp30(:,1:9);
        k_temp100=k_temp100(:,1:9);
        % 初始碳库
        % 用之前的平均最末尾的值，或者稳态值
        if tt ==1
            temp1 = wtchs30;
            temp2 = wtchs100;
        end              
        for i2 = Nt*(tt-1) + 2: Nt*tt+1
            outCmodel1 = Cmodel_vectorized(temp1, k_temp30);
            temp1 = outCmodel1 + temp1;
            outCmodel1 = Cmodel_vectorized(temp2, k_temp100);
            temp2 = outCmodel1 + temp2;
            % 如果当前年在目标列表里，就保存数据
            if ismember(i2/365 + 2020, target_years)
                Global_cpool_pre30(:,:,idd) = temp1 - wtchs30;
                Global_cpool_pre100(:,:,idd) = temp2 - wtchs100;
                idd = idd + 1 % 统一用 idd，修复你原来 dd 的bug！
            end 
        end
    end

% 20260524每年都有值
save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量全变year.mat', 'Global_cpool_pre30');
save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量全变year.mat', 'Global_cpool_pre100');

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量全变.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量全变.mat', 'Global_cpool_pre100');

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量只变NPP.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量只变NPP.mat', 'Global_cpool_pre100');

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量只变MAT.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量只变MAT.mat', 'Global_cpool_pre100');

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量只变MAP.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量只变MAP.mat', 'Global_cpool_pre100');

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量只变ET.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量只变ET.mat', 'Global_cpool_pre100');

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量只变AI.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量只变AI.mat', 'Global_cpool_pre100');


%% 计算SSP585在不同时间点的结果 
% 修改BaoHe3 的 3.2.2 部分，不区分SSP
% 3.2.2 将预测的未来的全球SSP的参数值代入模型 预测未来碳库  以及与历史的差值
% 修改第3.2部分，这部分之前没有输出中间的结果，重新再算一遍
% 修改用 ssp环境值拟合预测的历史环境值，再预测得到的只变一个环境的参数
% 使用ssp环境值拟合预测的历史环境值，再预测得到的全部用历史环境的参数  计算稳态值，作为未来预测的初始值
% 环境全变的SSP预测的参数 和原来相同  保持不变
indir='H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\';
indir2='H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数2\';  % 只变一个变量预测值
param_name = {'kbd_0_30_83', 'kbd_30_100_83', 'kl_0_30_83', 'kl_30_100_83', 'kp_0_30_83', 'kp_30_100_83', 'kq_0_30_83', 'kq_30_100_83', 'CUE'};
yt = {'2140','4160','6180','8100'};
% 读取SSP环境预测的历史环境 预测的历史参数值，用于计算稳态初始值
load([indir2,'Global_kbd_0_30_83_hist_all.mat']);
load([indir2,'Global_kl_0_30_83_hist_all.mat']);
load([indir2,'Global_kp_0_30_83_hist_all.mat']);
load([indir2,'Global_kq_0_30_83_hist_all.mat']);
load([indir2,'Global_CUE_hist_all.mat']);
load([indir2,'Global_kbd_30_100_83_hist_all.mat']);
load([indir2,'Global_kl_30_100_83_hist_all.mat']);
load([indir2,'Global_kp_30_100_83_hist_all.mat']);
load([indir2,'Global_kq_30_100_83_hist_all.mat']);

hh = 1837283;
param_SSP = single(zeros(hh,9,4)); % 只变一个变量预测值  9个变量 4个时间点
for ii =1:9
    for jj=1:4
%         files = [indir, 'Global_',param_name{ii},'_SSP585_',yt{jj},'.mat'];
%             files = [indir2, 'Global_',param_name{ii},'_SSP585_',yt{jj},'qu_MAT.mat']; % 只变一个变量预测值
%             files = [indir2, 'Global_',param_name{ii},'_SSP585_',yt{jj},'qu_MAP.mat']; % 只变一个变量预测值
%             files = [indir2, 'Global_',param_name{ii},'_SSP585_',yt{jj},'qu_ET.mat']; % 只变一个变量预测值
            files = [indir2, 'Global_',param_name{ii},'_SSP585_',yt{jj},'qu_AI.mat']; % 只变一个变量预测值
%         files = [indir2, 'Global_',param_name{ii},'_SSP585_',yt{jj},'qu_NPP.mat']; % 只变一个变量预测值
        temp = load(files);
        fieldNames = fieldnames(temp);
        originalVarName = fieldNames{1};
        temp2 = temp.(originalVarName);
        param_SSP(:,ii,jj)=temp2(:,3);
    end
end

NPP_SSP = single(zeros(1837283,4));
idh = 1;
for jj=1:4
	files = [indir, 'SSP585_',yt{jj},'的环境值.mat'];
	temp = load(files);
	fieldNames = fieldnames(temp);
	originalVarName = fieldNames{1};
	temp = temp.(originalVarName);
	NPP_SSP(:,idh)=temp(:,7);
	idh = idh+1;
end

% 只变一个变量预测值 取出对应的行=============
[~, ia, ~] = intersect(temp(:, 1:2), temp2(:, 1:2), 'rows');
NPP_SSP=NPP_SSP(ia,:);
temp=temp2;
% 历史NPP用拟合的值  
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\根据SSP线性拟合估计的历史环境值.mat')
NPP_HIST = final_result(ia,7);
% 计算稳态时的碳储量所需要的，就是历史NPP的碳输入
in0_30 = 3.474300E-03*NPP_HIST + 4.682669E-04;  % 任何时候都要计算这个
in30_100 = 1.342353E-03*NPP_HIST + 1.650560E-07;   % 任何时候都要计算这个
% 只变一个变量预测值 取出对应的行=============

% 准备其他碳库数据====
% 取出所有像元的经纬度
info = geotiffinfo('F:\黄土高原\参考蒸散量ET0数据\gaez-v5-data\temp.tif');
[lat,lon]=pix2latlon(info.RefMatrix,temp(:,1)+1,temp(:,2)+1); % 之前的索引从0开始的,在MATLAB中需要加1

% Qmax %  g C kg-1  用已知的栅格图，需要转为和参数一样的数据形式，一维
% 0.5度 由原来的数据直接进行单位转换，别的每没变
[Qmax0_30, R] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\Qmax0-30cm.tif');
[Qmax30_100, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\Qmax30-100cm.tif');
[row, col] = geographicToDiscrete(R, lat,lon);
Qmax0_30 = Qmax0_30(sub2ind(size(Qmax0_30), row, col)); %取出对应位置的值
Qmax30_100 = Qmax30_100(sub2ind(size(Qmax30_100), row, col)); %取出对应位置的值

% 常数值
pi=0.66;
kml = 0.001; % g C kg-1 day-1

% 运行时间 
years = 20;
Nt=years*365; % 年

k_jz=[kp_0_30_83_(:,3),kbd_0_30_83_(:,3),kq_0_30_83_(:,3),kl_0_30_83_(:,3),CUE_(:,3),in0_30,...
    repmat(pi, size(kp_0_30_83_, 1), 1),repmat(kml, size(kp_0_30_83_, 1), 1),Qmax0_30];
wtchs30 = wentai_vectorized(k_jz);
k_jz=[kp_30_100_83_(:,3),kbd_30_100_83_(:,3),kq_30_100_83_(:,3),kl_30_100_83_(:,3),CUE_(:,3),in30_100,...
    repmat(pi, size(kp_30_100_83_, 1), 1),repmat(kml, size(kp_30_100_83_, 1), 1),Qmax30_100];
wtchs100 = wentai_vectorized(k_jz);

Global_cpool_pre30 = single(zeros(hh,4,6));   % 点 4个库 6个时间点(2026 2039 2060 2078 2094 2100)   %最后一个时间点用于验证是否和之前算的一样
Global_cpool_pre100 = single(zeros(hh,4,6));   % 点 4个库 4个时间点(2026 2039 2060 2078 2094 2100)  %最后一个时间点用于验证是否和之前算的一样


% 定义你需要提取的目标年份
target_years = [2026, 2039, 2060, 2078, 2094, 2100];
idd=1; % 截取时间点计数器
% 每个时间段一个循环 21-40
    for tt = 1:4
        fprintf('year=%d\n',tt);
        kp30 = param_SSP(:,5,tt);
        kp100 = param_SSP(:,6,tt);
        kbd30 = param_SSP(:,1,tt);
        kbd100 = param_SSP(:,2,tt);
        kq30 = param_SSP(:,7,tt);
        kq100 = param_SSP(:,8,tt);
        kl30 = param_SSP(:,3,tt);
        kl100 = param_SSP(:,4,tt);
        cue = param_SSP(:,9,tt);
        %==========================
%         in0_30 = 3.474300E-03*NPP_SSP(:,tt) + 4.682669E-04;   % 全变和只变NPP的时候用这个数据
%         in30_100 = 1.342353E-03*NPP_SSP(:,tt) + 1.650560E-07;   % 全变和只变NPP的时候用这个数据
        %==========================
        % 修改后的参数
        k_temp30=[kp30,kbd30,kq30,kl30,cue,in0_30,repmat(pi, hh, 1),repmat(kml, hh, 1),Qmax0_30,wtchs30];
        k_temp100=[kp100,kbd100,kq100,kl100,cue,in30_100,repmat(pi, hh, 1),repmat(kml, hh, 1),Qmax30_100,wtchs100];  % 30-100
        nan_mask = any(isnan(k_temp30), 2);       % 标记含NaN的行
        neg_mask = any(k_temp30< 0, 2);          % 标记含负值的行
        target_rows = nan_mask | neg_mask; % 同时满足条件的行逻辑索引
        k_temp30(target_rows, :) = NaN;
        nan_mask = any(isnan(k_temp100), 2);       % 标记含NaN的行
        neg_mask = any(k_temp100< 0, 2);          % 标记含负值的行
        target_rows = nan_mask | neg_mask; % 同时满足条件的行逻辑索引
        k_temp100(target_rows, :) = NaN;
        
        k_temp30=k_temp30(:,1:9);
        k_temp100=k_temp100(:,1:9);
        % 初始碳库
        % 用之前的平均最末尾的值，或者稳态值
        if tt ==1
            temp1 = wtchs30;
            temp2 = wtchs100;
        end              
        for i2 = Nt*(tt-1) + 2: Nt*tt+1
            outCmodel1 = Cmodel_vectorized(temp1, k_temp30);
            temp1 = outCmodel1 + temp1;
            outCmodel1 = Cmodel_vectorized(temp2, k_temp100);
            temp2 = outCmodel1 + temp2;
            
            % 如果当前年在目标列表里，就保存数据
            if ismember(i2/365 + 2020, target_years)
                Global_cpool_pre30(:,:,idd) = temp1 - wtchs30;
                Global_cpool_pre100(:,:,idd) = temp2 - wtchs100;
                idd = idd + 1 % 统一用 idd，修复你原来 dd 的bug！
            end 
        end
        
    end

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量全变.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量全变.mat', 'Global_cpool_pre100');

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量只变NPP.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量只变NPP.mat', 'Global_cpool_pre100');

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量只变MAT.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量只变MAT.mat', 'Global_cpool_pre100');

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量只变MAP.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量只变MAP.mat', 'Global_cpool_pre100');

% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量只变ET.mat', 'Global_cpool_pre30');
% save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量只变ET.mat', 'Global_cpool_pre100');

save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量只变AI.mat', 'Global_cpool_pre30');
save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量只变AI.mat', 'Global_cpool_pre100');



%% 取出最后一年值做柱状图==================  hh,4,4,3 % 点 4个库 4个时间点 3个SSP 
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量全变.mat')
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量全变.mat')
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量只变MAT.mat')
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量只变MAT.mat')
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量只变MAP.mat')
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量只变MAP.mat')
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量只变AI.mat')
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量只变AI.mat')
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测0-30增量只变NPP.mat')
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数3\全球碳库SSP预测30-100增量只变NPP.mat')

load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量全变.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量全变.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量只变MAT.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量只变MAT.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量只变MAP.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量只变MAP.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量只变ET.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量只变ET.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量只变AI.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量只变AI.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测0-30增量只变NPP.mat')
load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\全球碳库SSP预测30-100增量只变NPP.mat')


Global_cpool_SSP30 = single(zeros(5*2,3,6));   %4个库+TOC库 平均值和标准差和95% 6个时间点
Global_cpool_SSP100 = single(zeros(5*2,3,6));   %4个库+TOC库 平均值和标准差和95% 6个时间点
for ssp=1:6  % 改为时间点一共6个
    temp = Global_cpool_pre30(:,:,ssp);  % 取最后一年
    clean_array = temp(all(~isinf(temp), 2),:);
    Global_cpool_SSP30(1:4,1,ssp) = mean(clean_array,'omitnan');
    Global_cpool_SSP30(1:4,2,ssp) = std(clean_array,1, 'omitnan');
    % 95%置信区间
    n = length(clean_array(~isnan(clean_array))); % 计算样本大小
    t_alpha_2_n_minus_1 = tinv(1-0.05/2, n-1); % 查找t分布的临界值，对于95%置信水平使用0.975（因为tinv需要双侧的概率）
    Global_cpool_SSP30(1:4,3,ssp) = t_alpha_2_n_minus_1 * (std(clean_array,1, 'omitnan') / sqrt(n));

    soc = sum(clean_array,2);
    Global_cpool_SSP30(5,1,ssp) = mean(soc,'omitnan');
    Global_cpool_SSP30(5,2,ssp) = std(soc,1, 'omitnan');
    % 95%置信区间
    n = length(clean_array(~isnan(clean_array))); % 计算样本大小
    t_alpha_2_n_minus_1 = tinv(1-0.05/2, n-1); % 查找t分布的临界值，对于95%置信水平使用0.975（因为tinv需要双侧的概率）
    Global_cpool_SSP30(5,3,ssp) = t_alpha_2_n_minus_1 * (std(soc,1, 'omitnan') / sqrt(n));

    % 去除异常值后的结果
    [out2,out3,out4] = time_yichang_2d(clean_array);
    [out5,out6,out7] = time_yichang_2d(soc);
    Global_cpool_SSP30(6:9,1,ssp) = out2;
    Global_cpool_SSP30(6:9,2,ssp) = out3;
    Global_cpool_SSP30(6:9,3,ssp) = out4;
    Global_cpool_SSP30(10,1,ssp) = out5;
    Global_cpool_SSP30(10,2,ssp) = out6;
    Global_cpool_SSP30(10,3,ssp) = out7;
    
    temp = Global_cpool_pre100(:,:,ssp);
    clean_array = temp(all(~isinf(temp), 2),:);
    Global_cpool_SSP100(1:4,1,ssp) = mean(clean_array,'omitnan');
    Global_cpool_SSP100(1:4,2,ssp) = std(clean_array,1, 'omitnan');
    % 95%置信区间
    n = length(clean_array(~isnan(clean_array))); % 计算样本大小
    t_alpha_2_n_minus_1 = tinv(1-0.05/2, n-1); % 查找t分布的临界值，对于95%置信水平使用0.975（因为tinv需要双侧的概率）
    Global_cpool_SSP100(1:4,3,ssp) = t_alpha_2_n_minus_1 * (std(clean_array,1, 'omitnan') / sqrt(n));

    soc = sum(clean_array,2);
    Global_cpool_SSP100(5,1,ssp) = mean(soc,'omitnan');
    Global_cpool_SSP100(5,2,ssp) = std(soc,1, 'omitnan');
    % 95%置信区间
    n = length(clean_array(~isnan(clean_array))); % 计算样本大小
    t_alpha_2_n_minus_1 = tinv(1-0.05/2, n-1); % 查找t分布的临界值，对于95%置信水平使用0.975（因为tinv需要双侧的概率）
    Global_cpool_SSP100(5,3,ssp) = t_alpha_2_n_minus_1 * (std(soc,1, 'omitnan') / sqrt(n));

    % 去除异常值后的结果
    [out2,out3,out4] = time_yichang_2d(clean_array);
    [out5,out6,out7] = time_yichang_2d(soc);
    Global_cpool_SSP100(6:9,1,ssp) = out2;
    Global_cpool_SSP100(6:9,2,ssp) = out3;
    Global_cpool_SSP100(6:9,3,ssp) = out4;
    Global_cpool_SSP100(10,1,ssp) = out5;
    Global_cpool_SSP100(10,2,ssp) = out6;
    Global_cpool_SSP100(10,3,ssp) = out7;
end

%4个库+TOC库 平均值和标准差和95% 6个时间点
pjbi = 3; % 指定取标准差还是95%
zs = 0; % 下面是未去除异常值的结果
temp = single(zeros(30,8));
temp([2:5,7:10,12:15,17:20,22:25,27:30],1) = reshape(squeeze(Global_cpool_SSP30(1+zs:4+zs,1,:)), [], 1); % 分库
temp([2:5,7:10,12:15,17:20,22:25,27:30],2) = reshape(squeeze(Global_cpool_SSP30(1+zs:4+zs,pjbi,:)), [], 1); % 分库
temp([1,6,11,16,21,26],1) = reshape(squeeze(Global_cpool_SSP30(5+zs,1,:)), [], 1); % 分库
temp([1,6,11,16,21,26],2) = reshape(squeeze(Global_cpool_SSP30(5+zs,pjbi,:)), [], 1); % 分库
temp([2:5,7:10,12:15,17:20,22:25,27:30],3) = reshape(squeeze(Global_cpool_SSP100(1+zs:4+zs,1,:)), [], 1); % 分库
temp([2:5,7:10,12:15,17:20,22:25,27:30],4) = reshape(squeeze(Global_cpool_SSP100(1+zs:4+zs,pjbi,:)), [], 1); % 分库
temp([1,6,11,16,21,26],3) = reshape(squeeze(Global_cpool_SSP100(5+zs,1,:)), [], 1); % 分库
temp([1,6,11,16,21,26],4) = reshape(squeeze(Global_cpool_SSP100(5+zs,pjbi,:)), [], 1); % 分库
zs = 5; % 下面是去除异常值的结果
temp([2:5,7:10,12:15,17:20,22:25,27:30],5) = reshape(squeeze(Global_cpool_SSP30(1+zs:4+zs,1,:)), [], 1); % 分库
temp([2:5,7:10,12:15,17:20,22:25,27:30],6) = reshape(squeeze(Global_cpool_SSP30(1+zs:4+zs,pjbi,:)), [], 1); % 分库
temp([1,6,11,16,21,26],5) = reshape(squeeze(Global_cpool_SSP30(5+zs,1,:)), [], 1); % 分库
temp([1,6,11,16,21,26],6) = reshape(squeeze(Global_cpool_SSP30(5+zs,pjbi,:)), [], 1); % 分库
temp([2:5,7:10,12:15,17:20,22:25,27:30],7) = reshape(squeeze(Global_cpool_SSP100(1+zs:4+zs,1,:)), [], 1); % 分库
temp([2:5,7:10,12:15,17:20,22:25,27:30],8) = reshape(squeeze(Global_cpool_SSP100(1+zs:4+zs,pjbi,:)), [], 1); % 分库
temp([1,6,11,16,21,26],7) = reshape(squeeze(Global_cpool_SSP100(5+zs,1,:)), [], 1); % 分库
temp([1,6,11,16,21,26],8) = reshape(squeeze(Global_cpool_SSP100(5+zs,pjbi,:)), [], 1); % 分库


% 直接读取结果数据，做掩膜，再统计
[TOC301_5, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_1.5C增量TOC.tif');
[EAS, ~] = geotiffread('H:\数据库\ESA_CCI_LandCoverData\重采样.tif');  % 掩膜掉10 11 12 20 30 190 200 201 202 210 220
[KOP, ~] = geotiffread('H:\数据库\气候类型数据\koppen_geiger_tif\1991_2020\重采样.tif');  % 掩膜掉29
% 条件1：TOC301_5为NaN
cond1 = isnan(TOC301_5);

% 条件2：EAS属于指定类别
eas_list = [10,11,12,20,30,190,200,201,202,210,220];
cond2 = ismember(EAS, eas_list);

% 条件3：KOP等于29
cond3 = (KOP == 29);

% 三者并集：任意条件成立则掩膜位置为0
mask_union = cond1 | cond2 | cond3;

% 生成最终掩膜：mask=0为剔除，mask=1为保留
mask = ones(size(TOC301_5), 'uint8');
mask(mask_union) = 0;

% [TOC301_5, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_1.5C增量TOC.tif');
% [TOC302, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_2.0C增量TOC.tif');
% [TOC303, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_3.0C增量TOC.tif');
% [TOC304, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_4.0C增量TOC.tif');
% [TOC305, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_5.0C增量TOC.tif');
% 换为四个碳库
[TOC301_5, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_1.5C增量.tif');
[TOC302, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_2.0C增量.tif');
[TOC303, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_3.0C增量.tif');
[TOC304, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_4.0C增量.tif');
[TOC305, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\增量碳库预测栅格4\全球碳库SSP预测30_5.0C增量.tif');
TOC301_5 = TOC301_5(:,:,4);
TOC302 = TOC302(:,:,4);
TOC303 = TOC303(:,:,4);
TOC304 = TOC304(:,:,4);
TOC305 = TOC305(:,:,4);
TOC301_5 = TOC301_5 .* double(mask);
TOC302 = TOC302 .* double(mask);
TOC303 = TOC303 .* double(mask);
TOC304 = TOC304 .* double(mask);
TOC305 = TOC305 .* double(mask);
TOC301_5(TOC301_5 == 0) = NaN;
TOC302(TOC302 == 0) = NaN;
TOC303(TOC303 == 0) = NaN;
TOC304(TOC304 == 0) = NaN;
TOC305(TOC305 == 0) = NaN;

% 将每个栅格展成一维列向量
v1 = TOC301_5(~isnan(TOC301_5));
v2 = TOC302(~isnan(TOC302));
v3 = TOC303(~isnan(TOC303));
v4 = TOC304(~isnan(TOC304));
v5 = TOC305(~isnan(TOC305));

mu0 = mean(v1);
sigma0 = std(v1);
lower_thr = mu0 - 3*sigma0;
upper_thr = mu0 + 3*sigma0;
v1_clean = v1(v1 >= lower_thr & v1 <= upper_thr);
mu0 = mean(v2);
sigma0 = std(v2);
lower_thr = mu0 - 3*sigma0;
upper_thr = mu0 + 3*sigma0;
v2_clean = v2(v2 >= lower_thr & v2 <= upper_thr);
mu0 = mean(v3);
sigma0 = std(v3);
lower_thr = mu0 - 3*sigma0;
upper_thr = mu0 + 3*sigma0;
v3_clean = v3(v3 >= lower_thr & v3 <= upper_thr);
mu0 = mean(v4);
sigma0 = std(v4);
lower_thr = mu0 - 3*sigma0;
upper_thr = mu0 + 3*sigma0;
v4_clean = v4(v4 >= lower_thr & v4 <= upper_thr);
mu0 = mean(v5);
sigma0 = std(v5);
lower_thr = mu0 - 3*sigma0;
upper_thr = mu0 + 3*sigma0;
v5_clean = v5(v5 >= lower_thr & v5 <= upper_thr);

% 拼接全部数据
data_all = [v1_clean; v2_clean; v3_clean; v4_clean; v5_clean];

% 生成分组标签，数字1?5对应5个箱子
g1 = ones(length(v1_clean),1);
g2 = 2*ones(length(v2_clean),1);
g3 = 3*ones(length(v3_clean),1);
g4 = 4*ones(length(v4_clean),1);
g5 = 5*ones(length(v5_clean),1);
group = [g1;g2;g3;g4;g5];

% MATLAB端：保存mat
save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\maoc_groups.mat',...
    'data_all','group', '-v6');
save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\mbc_groups.mat',...
    'data_all','group', '-v6');
save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\doc_groups.mat',...
    'data_all','group', '-v6');
save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\poc_groups.mat',...
    'data_all','group', '-v6');
save('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\未来情景的预测值\全球只变单一变量预测的参数4\toc_groups.mat',...
    'data_all','group', '-v6');


figure('Color','w');
bp = boxplot(data_all, group, ...
    'Labels',{'TOC301','TOC302','TOC303','TOC304','TOC305'},...
    'Whisker',1.5,...
    'Symbol','');   %隐藏离群点

xlabel('Dataset');
ylabel('TOC value');
title('Boxplot of five TOC datasets');
grid on;
hold on;

% 替换原来mean()，用10%截尾均值，去掉上下各10%极端值
trim_percent = 15;   %上下各剔除10%，可改为5

mu_trim(1) = trimmean(v1_clean, trim_percent);
mu_trim(2) = trimmean(v2_clean, trim_percent);
mu_trim(3) = trimmean(v3_clean, trim_percent);
mu_trim(4) = trimmean(v4_clean, trim_percent);
mu_trim(5) = trimmean(v5_clean, trim_percent);

% IQR四分位法 对5组分别剔除异常
v1_noout = iqr_remove_outlier(v1_clean);
v2_noout = iqr_remove_outlier(v2_clean);
v3_noout = iqr_remove_outlier(v3_clean);
v4_noout = iqr_remove_outlier(v4_clean);
v5_noout = iqr_remove_outlier(v5_clean);

% 在剔除异常后的数据集上算普通均值
mu_noout(1)=mean(v1_noout);
mu_noout(2)=mean(v2_noout);
mu_noout(3)=mean(v3_noout);
mu_noout(4)=mean(v4_noout);
mu_noout(5)=mean(v5_noout);

xpos = 1:5;   %5个箱子x坐标位置

% 叠加均值，黑色菱形标记，Marker可以改：'o'圆形, '^'三角
plot(xpos, mu_trim, 'kd','MarkerSize',7,'MarkerFaceColor','k','DisplayName','Mean');

% 箱体美化填充
h = findobj(gca,'Tag','Box');
for i=1:length(h)
    patch(get(h(i),'XData'),get(h(i),'YData'),[0.6 0.75 0.9],'FaceAlpha',0.6);
end

legend('Location','best');
set(gca,'FontName','Arial','FontSize',11);
hold off;

% 绘制箱线图
figure('Color','w');
boxplot(data_all, group, ...
    'Labels',{'TOC301','TOC302','TOC303','TOC304','TOC305'},...
    'Whisker',1.5,...
    'Symbol','');   % =========关键：隐藏所有离群值标记=========

xlabel('Dataset');
ylabel('TOC value');
title('Boxplot of five TOC datasets');
grid on;


% 使用中位数方法检测并删除离群值
vec_all = TOC305(~isnan(TOC305));   % 提取全部非NaN，转为一维向量
[out5,out6,out7] = time_yichang_2d(vec_all);

% 按照下面的方法统计 一直是减少的，而且没有明显趋势，
% 三倍标准差 去除的异常值比较不错但是结果不好
vec_all = TOC301_5(~isnan(TOC301_5));   % 提取全部非NaN，转为一维向量
if length(vec_all) <2
    mean_val=NaN;ci_low=NaN;ci_high=NaN;n_after=0;
else
    mu0 = mean(vec_all);
    sigma0 = std(vec_all);
    lower_thr = mu0 - 3*sigma0;
    upper_thr = mu0 + 3*sigma0;

    vec_clean = vec_all(vec_all >= lower_thr & vec_all <= upper_thr);
    n_all = length(vec_all);
    n_after = length(vec_clean);

    if n_after<2
        mean_val=NaN;ci_low=NaN;ci_high=NaN;
    else
        mean_val = mean(vec_clean);
        std_val = std(vec_clean);
        se = std_val / sqrt(n_after);
        t_crit = tinv(0.975, n_after-1);
        ci_low = mean_val - t_crit*se;
        ci_high = mean_val + t_crit*se;
    end
end

fprintf('原始n=%d，去异常后n=%d\n',n_all,n_after);
fprintf('3?sigma阈值 [%.4f,%.4f]\n',lower_thr,upper_thr);
fprintf('均值=%.4f，95%%CI [%.4f,%.4f]\n',mean_val,ci_low,ci_high);


%% 数据验证
temp = readtable('H:\数据库\基于青年基金收集的数据\Manuscript\SOCTimeSeries-main\SOC_withNewData_for_analysis_2025Jul21.csv');
soc=temp(:,4:7);
soc = table2array(soc);
soc = [(1:length(soc))', soc];  
% 1. 按经纬度分组，给相同位置分配相同站点ID
% 提取经纬度配对
lonlat = soc(:, [3, 2]);  % [lon, lat] 配对，用于分组
% 对经纬度去重，获取唯一站点
[~, ~, site_id] = unique(lonlat, 'rows');  
% 2. 拼接新矩阵：新增站点ID为第一列
soc = [soc(:,1),site_id, soc(:,2:end)];  
clear temp

%% 将所有的数据合并 20260421
npp2000_2025 = 'H:\数据库\基于青年基金收集的数据\Manuscript\SOCTimeSeries-main\MODIS_NPP_WideTable.mat';
npp1982_2000 = 'H:\数据库\基于青年基金收集的数据\Manuscript\SOCTimeSeries-main\GLASS_NPP_site_extract.mat';
npp1943_2022 = 'H:\碳转移效率全球研究\data\Processed_NPP_Annual\NPP_41276points_4indices_80years.mat';
load(npp2000_2025);
npp2000_2025 = output;
load(npp1982_2000);
npp1982_2000 = npp_site_array;
load(npp1943_2022);
npp1943_2022 = output;

output = zeros(41276, 3 + 83);
output(:,:) = NaN;  % 全部先设为空
output(:,1) = (1:41276);
output(:,2) = lat_sites;
output(:,3) = lon_sites;
% 矫正数据
output(:,4:42)  = npp1943_2022(:,1,1:39);    % 1943-2022 (39年)
output(:,43:61) = npp1982_2000(:,1:19);      % 1982-2000 (19年)
output(:,62:86) = npp2000_2025(:,5:29);     % 2001-2025 (25年) → 基准

output2=output(:,4:end);
% 统一无效值：0 → NaN
output2(output2 == 0) = NaN;
output2(output2 == 65.5350) = NaN;
output2(output2 < 0) = NaN;
output2 = [output(:,1:3), output2];

load('H:\数据库\基于青年基金收集的数据\Manuscript\SOCTimeSeries-main\Jia_site_param.mat');
% 求出 A、B、C 中 任意一行含 NaN 的行号 并集
nan_rows = unique( [ find(any(isnan(output2),2)) ; find(any(isnan(Jia_site_param),2));] );
% 删除并集中的所有行
soc(nan_rows,:) = [];  % SOC

% 处理SOC数据
% ===================== 找出每个站点有多少年 =====================
siteIDs = unique(soc(:,2));        % 所有站点ID
nSites = length(siteIDs);

% 统计每个站点的年份数量
yearCounts = zeros(nSites,1);

for i = 1:nSites
    idx = soc(:,2) == siteIDs(i);
    yearCounts(i) = length(unique(soc(idx,5)));  % 去重统计年份
end

% ===================== 3. 找出“只有1个年份”的站点 =====================
badSites = siteIDs(yearCounts == 1);   % 只含1年的站点
goodSites = siteIDs(yearCounts >= 2);  % 多年份站点

fprintf('总站点数：%d\n', nSites);
fprintf('只含1年的站点数：%d\n', length(badSites));
fprintf('保留多年份站点数：%d\n', length(goodSites));

% ===================== 4. 过滤数据（保留多年份站点） =====================
% 标记哪些行属于多年份站点
keepRows = ismember(soc(:,2), goodSites);

% 过滤后的数据
soc_filtered = soc(keepRows,:);

% ===================== 5. 结果统计 =====================
fprintf('原始数据行数：%d\n', size(soc,1));
fprintf('过滤后行数：%d\n', size(soc_filtered,1));
fprintf('删除行数：%d\n', size(soc,1) - size(soc_filtered,1));


% 同一站点同一年份有多个观测值的计算平均值标准差
% 提取分组键：【站点ID, 纬度, 经度, 年份】四列完全一致才为同一组
group_key = soc_filtered(:, [2,3,4,5]);

% 获取所有唯一分组，以及每行对应的分组编号
[unique_group, ~, group_idx] = unique(group_key, 'rows');
n_unique_group = size(unique_group, 1);

% 预分配最终输出矩阵：7列 [全局ID, 站点ID, lat, lon, year, 均值, 标准差]
soc_final = zeros(n_unique_group, 7);

% 遍历每一个唯一分组，依次填充
for g = 1:n_unique_group
    % 找到当前分组对应的所有原始行
    rows_in_group = (group_idx == g);
    
    % 1. 第1列：保留该组【第一个原始全局ID】
    soc_final(g, 1) = soc_filtered(find(rows_in_group, 1), 1);
    
    % 2. 第2-5列：保留该组唯一值（站点ID、lat、lon、year，组内全重复）
    soc_final(g, 2:5) = unique_group(g, :);
    
    % 3. 提取该组所有SOC原始值
    soc_vals = soc_filtered(rows_in_group, 6);
    
    % 4. 第6列：均值（自动忽略NaN）
    soc_final(g, 6) = mean(soc_vals, 'omitnan');
    
    % 5. 第7列：标准差（自动忽略NaN）
    soc_final(g, 7) = std(soc_vals, 'omitnan');
end

soc_final = sortrows(soc_final, 1);
% 结果信息输出
fprintf('===== 聚合去重完成 =====\n');
fprintf('原始数据行数：%d\n', size(soc_filtered, 1));
fprintf('去重聚合后行数：%d\n', size(soc_final, 1));
fprintf('最终列顺序：\n');
fprintf('1:全局ID(组首值) | 2:站点ID | 3:纬度 | 4:经度 | 5:年份 | 6:SOC均值 | 7:SOC标准差\n');


% 取出对应位置的NPP和参数值
% 矫正NPP
output2=output2(soc_final(:,1),:);  % NPP
npp1943 = squeeze(npp1943_2022(soc_final(:,1),1,40:44)); % 1982年
NPP=output2;
NPP(:,4:42) = output2(:,4:42).*mean(output2(:,43:47),2)./mean(npp1943,2);
% 参数
Jia_site_param = Jia_site_param(soc_final(:,1),:); 

% 根据已知的I值与NPP的关系计算新的I值，用线性拟合
aa=[mean(NPP(:,63:84),2),Jia_site_param(:,12)];
aa(aa(:,2)>0.01,:)=[];
p = polyfit(aa(:,1),aa(:,2),1);
NPP2=polyval(p,NPP(:,4:end));  % 估计的I值


% 模拟
% ===================== 0. 全局基础参数=====================
pi  = 0.66;
kml = 0.001;    % g C kg-1 day-1
total_year = 83;% 完整时间轴：1943~2025 共83年
BD30 = mean(Jia_site_param(:,[24,24,25]),2);   % 土壤容重 g cm-3

% ===================== 1. 提取所有唯一站点 & 绑定站点全部信息 =====================
% 输入：soc_final = [全局ID, 站点ID, lat, lon, year, SOC均值, SOC标准差]
all_siteID = unique(soc_final(:,2));   % 全部不重复站点
site_num = length(all_siteID);         % 总站点数量

% 预分配三维模拟数组：【站点数 × 4个碳库 × 83年】
cpoolarr30 = single(zeros(site_num, 4, total_year)); 
site_id = zeros(site_num,1); % 站点ID
BD30_site = zeros(site_num,1); % 站点BD的唯一值
% ===================== 2. 逐站点绑定：最早年份初始值 + 逐年观测SOC/NPP =====================
year2idx = @(year) year - 1942;  

for site_i = 1:site_num
    % 当前遍历站点ID
    now_site = all_siteID(site_i);
    site_id(site_i) = now_site;
    % 提取该站点所有观测行
    site_data = soc_final(soc_final(:,2)==now_site, :);
    aa = BD30(soc_final(:,2)==now_site, :);
    BD30_site(site_i) = aa(1);
    % ---------- 2.1 提取该站点最早观测年份 & 对应SOC初始值 ----------
    site_year_all = site_data(:,5);
    site_earliest_year = min(site_year_all);          % 站点最早观测年
    init_idx = year2idx(site_earliest_year);          % 对应时间维序号
    init_SOC_mean = site_data(site_year_all==site_earliest_year, 6); % 初始SOC均值
    
    % ---------- 2.2 沿用原始4碳库初始分配规则 ----------
    % 初始库：Jia_site_param 对应站点的4个碳库参数，绑定站点ID匹配
    site_param = Jia_site_param(soc_final(:,2)==now_site, :);
    % ===================== 这里开始：按比例分配初始值 =====================
    original_pools = site_param(1, [16,18,20,22]);
    total_original = sum(original_pools);
    if total_original > 0
        pool_ratios = original_pools / total_original;
    else
        pool_ratios = [0.25, 0.25, 0.25, 0.25];
    end
    
    % 分配新初始库
    % 转换：Mg C ha?? （30cm土层深度、容重加权）
    init_SOC_mean = init_SOC_mean ./ (BD30_site(site_i) * 30 * 0.1);
    new_initial_pools = init_SOC_mean * pool_ratios;
    cpoolarr30(site_i, :, init_idx) = new_initial_pools;
    
    % ---------- 2.3 预存该站点所有年份的观测SOC、对应NPP输入 ----------
    % 建立【年份→SOC均值、年份→NPP】映射表
    site_year_list = site_data(:,5);
    site_SOC_obs = site_data(:,6);
    % 匹配对应年份的站点NPP（NPP矩阵站点、年份完全一一对应）
    NPP3 = NPP2(soc_final(:,2)==now_site,:);
    
    % 时间循环：从初始年份的下一年，一直迭代到最终2025年（第83年）
    for t = init_idx+1: max(year2idx(site_year_list))
        % 上一年碳库状态
        cpool_last = cpoolarr30(site_i,:,t-1);
    
        % 逐年NPP驱动输入
        NPP_now = NPP3(1, t);   % 每个站点第t年对应的真实NPP值
    
        % 完全沿用原始模型函数输入格式
        model_input = [site_param(1,[7,3,9,5,11]),...
            NPP_now,...   % 替换为逐年真实NPP驱动！！！
            pi,...
            kml,...
            site_param(1,14)];
    
        % 模型计算年新增碳输入
        outCmodel1 = Cmodel_vectorized(cpool_last, model_input);
    
        % 碳库累加迭代（和原版迭代逻辑完全一致）
        cpoolarr30(site_i,:,t) = cpool_last + outCmodel1;
    end
end

% ===================== 4. 单位换算（100%沿用你原始公式）=====================
% 转换：Mg C ha?? （30cm土层深度、容重加权）
cpoolarr30 = cpoolarr30 .* BD30_site * 30 * 0.1;

% ===================== 5. 输出总碳库（4库求和，二维矩阵：站点 × 年份）=====================
moni_cpool = squeeze(sum(cpoolarr30, 2));  

% 有一行异常值去掉
moni_cpool(869,56:end) = 0;


% 自动提取所有站点观测年份的模拟SOC
% 这个包含初始值
moni_obs_match = [];
for site_i = 1:site_num
    now_site = all_siteID(site_i);
    site_data = soc_final(soc_final(:,2)==now_site, :);
    site_year = site_data(:,5);
    y_idx = year2idx(site_year);
    % 该站点所有观测年的模拟SOC
    moni_val = moni_cpool(site_i, y_idx);
    obs_val = site_data(:,6);
    moni_obs_match = [moni_obs_match; repmat(now_site,length(obs_val),1), site_year, obs_val, moni_val'];
end

% 这个不包含初始值
moni_obs_match2 = [];
for site_i = 1:site_num
    now_site = all_siteID(site_i);
    site_data = soc_final(soc_final(:,2)==now_site, :);
    
    site_year_all = site_data(:,5);
    site_earliest_year = min(site_year_all);  % 初始年份（要剔除）
    
    % ===================== 剔除初始年份 =====================
    keep_idx = site_year_all ~= site_earliest_year;  % 只保留非初始年
    site_data_keep = site_data(keep_idx, :);
    
    if isempty(site_data_keep)
        continue; % 没有可保留数据就跳过
    end
    
    site_year = site_data_keep(:,5);
    obs_val  = site_data_keep(:,6);
    y_idx = year2idx(site_year);
    
    % 该站点 非初始年 的模拟SOC
    moni_val = moni_cpool(site_i, y_idx);
    
    % 拼接结果
    moni_obs_match2 = [moni_obs_match2; 
        repmat(now_site,length(obs_val),1), ...
        site_year, ...
        obs_val, ...
        moni_val'];
end


figure('Color','w'); hold on; grid on;
all_site = unique(moni_obs_match(:,1));  % 提取所有唯一站点ID
% 颜色循环（自动给每个站点分配不同颜色，区分清晰）
color_list = lines(length(all_siteID));
for i = 1:100    %length(all_site)
    sid = all_site(i);
    % 提取当前站点所有数据
    site_data = moni_obs_match(moni_obs_match(:,1)==sid, :);
    year_obs = site_data(:,2);    % 实测年份
    obs = site_data(:,3);         % 实测SOC（散点）
    y_idx = year2idx(year_obs);
    
    % ========== 提取该站点完整全时序模拟线（1943-2025） ==========
    site_full_sim = moni_cpool(i, min(y_idx):max(y_idx));  % 该站点完整83年模拟序列
    full_year = (min(year_obs):max(year_obs));             % 全域时间轴
    
    % ========== 绘图 ==========
    % 1. 模拟完整趋势线（整条时间序列红线）
    plot(full_year, site_full_sim, '-', 'Color', color_list(i,:), 'LineWidth',1.5, 'HandleVisibility','off');
    % 2. 实测观测散点（空心点，醒目区分）
    scatter(year_obs, obs, 40, color_list(i,:), 'filled', 'DisplayName', sprintf('站点%d',sid));
end

% 图表美化
xlabel('年份','FontSize',12);
ylabel('SOC 含量 (Mg C ha^{-1})','FontSize',12);
title('各站点SOC实测值(散点)与模型模拟全时序趋势(曲线)','FontSize',13);
legend('Location','best');
xlim([1940 2030]);  % 年份轴范围适配
hold off;

% 保存到当前文件夹
save('H:\数据库\基于青年基金收集的数据\Manuscript\SOCTimeSeries-main\包含初始值的对应时间模拟.mat', 'moni_obs_match', '-v6');
save('H:\数据库\基于青年基金收集的数据\Manuscript\SOCTimeSeries-main\不包含初始值的对应时间模拟.mat', 'moni_obs_match2', '-v6');
save('H:\数据库\基于青年基金收集的数据\Manuscript\SOCTimeSeries-main\每个站点在观测的时间内的模拟结果.mat', 'moni_cpool', '-v6');
save('H:\数据库\基于青年基金收集的数据\Manuscript\SOCTimeSeries-main\站点ID.mat', 'site_id', '-v6');


%% 20260411====================
% 1，首先提取对应位置的历史参数值
% 以0.008333度分辨率，0-30 30-100间隔分层的结果
% load('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\随机森林相对重要性\随机选择点的参数和环境值_已去除异常值2.mat');
latlon= soc(:,1:2);
param2 = zeros(size(latlon(:,1),1),36);
param2(:,1:2) = latlon;
% 读取数据
[kbd_0_30, R] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\0-30cm与30-100cm间隔分层0.008333度\kbd_0_30.tif');
[row, col] = geographicToDiscrete(R, latlon(:,1), latlon(:,2));  % 转换为栅格行列号
param2(:,3) = kbd_0_30(sub2ind(size(kbd_0_30), row, col)); %取出对应位置的值
clear kbd_0_30
[kbd_30_100, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\0-30cm与30-100cm间隔分层0.008333度\kbd_30_100.tif');
param2(:,4) = kbd_30_100(sub2ind(size(kbd_30_100), row, col)); %取出对应位置的值
clear kbd_30_100
[kl_0_30, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\0-30cm与30-100cm间隔分层0.008333度\kl_0_30.tif');
param2(:,5) = kl_0_30(sub2ind(size(kl_0_30), row, col)); %取出对应位置的值
clear kl_0_30
[kl_30_100, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\0-30cm与30-100cm间隔分层0.008333度\kl_30_100.tif');
param2(:,6) = kl_30_100(sub2ind(size(kl_30_100), row, col)); %取出对应位置的值
clear kl_30_100
[kp_0_30, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\0-30cm与30-100cm间隔分层0.008333度\kp_0_30.tif');
param2(:,7) = kp_0_30(sub2ind(size(kp_0_30), row, col)); %取出对应位置的值
clear kp_0_30
[kp_30_100, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\0-30cm与30-100cm间隔分层0.008333度\kp_30_100.tif');
param2(:,8) = kp_30_100(sub2ind(size(kp_30_100), row, col)); %取出对应位置的值
clear kp_30_100
[kq_0_30, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\0-30cm与30-100cm间隔分层0.008333度\kq_0_30.tif');
param2(:,9) = kq_0_30(sub2ind(size(kq_0_30), row, col)); %取出对应位置的值
clear kq_0_30
[kq_30_100, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\0-30cm与30-100cm间隔分层0.008333度\kq_30_100.tif');
param2(:,10) = kq_30_100(sub2ind(size(kq_30_100), row, col)); %取出对应位置的值
clear kq_30_100


% kbd 计算为kbd * MBC
% 只按0-30  30-100cm分层 单位 mgC kg-1Soil;
[mbc30_100, R1] = geotiffread('j:\组里数据备份-2023.12.18\global pattern\MBC_30-100cm.tif');    
[mbc0_30, R2] = geotiffread( 'j:\组里数据备份-2023.12.18\global pattern\MBC30cm.tif');
mask2 = (mbc30_100 <= 0) | (mbc0_30 <= 0);
mbc30_100(mask2) = NaN;
mbc0_30(mask2) = NaN;
mbc30_100 = mbc30_100/1000; % gC kg-1Soil;
mbc0_30 = mbc0_30/1000; % gC kg-1Soil;
clear mask2
[row, col] = geographicToDiscrete(R1, latlon(:,1), latlon(:,2));  % 转换为栅格行列号
param2_mbc30_100 = mbc30_100(sub2ind(size(mbc30_100), row, col)); %取出对应位置的值 
[row, col] = geographicToDiscrete(R2, latlon(:,1), latlon(:,2));  % 转换为栅格行列号
param2_mbc0_30 = mbc0_30(sub2ind(size(mbc0_30), row, col)); %取出对应位置的值 
param2(:,3)=param2(:,3).*param2_mbc0_30;
param2(:,4)=param2(:,4).*param2_mbc30_100;
clear mbc0_30 mbc30_100 R1 R2

% 加入CUE变量
% CUE 这个只用0-30cm单层的, 和MBC的分辨率一样
[CUE, R] = geotiffread('J:\组里数据备份-2023.12.18\global pattern\CUE_y.tif');
CUE(CUE<=0) = NaN;
[row, col] = geographicToDiscrete(R, latlon(:,1), latlon(:,2));  % 转换为栅格行列号
param2_cue = CUE(sub2ind(size(CUE), row, col)); %取出对应位置的值
param2(:,11) = param2_cue;
clear param2_cue CUE


% 土壤性质数据
% Sand Silt Clay texture BD SOC_gkg TN CN PH CEC
filename = 'F:\Google下载\wise_30sec_v1\wise_30sec_v1\导出实验\wise30sec_fin.tif';
[bdgis, R] = geotiffread(filename);
[row, col] = geographicToDiscrete(R, latlon(:,1), latlon(:,2));  % 转换为栅格行列号
bdgis_value = bdgis(sub2ind(size(bdgis), row, col)); %取出对应位置的值
clear bdgis 
[wise30,NEWSUID,~] =xlsread('F:\Google下载\wise_30sec_v1\wise_30sec_v1\导出实验\wise30sec_fin.xlsx','wise30sec_fin');
NEWSUID = string(NEWSUID);
CODE = NEWSUID(2:end,5);
CODE(CODE == 'A') =1;
CODE(CODE == 'B') =2;
CODE(CODE == 'C') =3;
CODE(CODE == 'D') =4;
CODE(CODE == 'E') =5;

[wD1,NEWSUID,~] =xlsread('F:\Google下载\wise_30sec_v1\wise_30sec_v1\导出实验\HW30s_wD1_to_excel.xls','HW30s_wD1_to_excel');
NEWSUID = string(NEWSUID);
pscl_wD1 = NEWSUID(2:end,23);
% 多重替换 # 将类别:C  M  F  V  Z -  设置为：1  2  3  4  5 0
pscl_wD1 = strrep(pscl_wD1, 'M', '2');
pscl_wD1 = strrep(pscl_wD1, 'F', '3'); 
pscl_wD1 = strrep(pscl_wD1, 'C', '1');
pscl_wD1 = strrep(pscl_wD1, '-', '0');
pscl_wD1 = strrep(pscl_wD1, 'V', '4');
pscl_wD1 = strrep(pscl_wD1, 'Z', '5');
% 转换为数值数组
pscl_wD1 = str2double(pscl_wD1);
NEWSUID_wD1 = NEWSUID(2:end,1);

[wD2,NEWSUID,~] =xlsread('F:\Google下载\wise_30sec_v1\wise_30sec_v1\导出实验\HW30s_wD2_to_excel.xls','HW30s_wD2_to_excel');
NEWSUID = string(NEWSUID);
pscl_wD2 = NEWSUID(2:end,23);
% 多重替换
pscl_wD2 = strrep(pscl_wD2, 'M', '2');
pscl_wD2 = strrep(pscl_wD2, 'F', '3'); 
pscl_wD2 = strrep(pscl_wD2, 'C', '1');
pscl_wD2 = strrep(pscl_wD2, '-', '0');
pscl_wD2 = strrep(pscl_wD2, 'V', '4');
% 转换为数值数组
pscl_wD2 = str2double(pscl_wD2);
NEWSUID_wD2 = NEWSUID(2:end,1);

[wD3,NEWSUID,~] =xlsread('F:\Google下载\wise_30sec_v1\wise_30sec_v1\导出实验\HW30s_wD3_to_excel.xls','HW30s_wD3_to_excel');
NEWSUID = string(NEWSUID);
pscl_wD3 = NEWSUID(2:end,23);
% 多重替换
pscl_wD3 = strrep(pscl_wD3, 'M', '2');
pscl_wD3 = strrep(pscl_wD3, 'F', '3'); 
pscl_wD3 = strrep(pscl_wD3, 'C', '1');
pscl_wD3 = strrep(pscl_wD3, '-', '0');
pscl_wD3 = strrep(pscl_wD3, 'V', '4');
% 转换为数值数组
pscl_wD3 = str2double(pscl_wD3);
NEWSUID_wD3 = NEWSUID(2:end,1);

[wD4,NEWSUID,~] =xlsread('F:\Google下载\wise_30sec_v1\wise_30sec_v1\导出实验\HW30s_wD4_to_excel.xls','HW30s_wD4_to_excel');
NEWSUID = string(NEWSUID);
pscl_wD4 = NEWSUID(2:end,23);
% 多重替换
pscl_wD4 = strrep(pscl_wD4, 'M', '2');
pscl_wD4 = strrep(pscl_wD4, 'F', '3'); 
pscl_wD4 = strrep(pscl_wD4, 'C', '1');
pscl_wD4 = strrep(pscl_wD4, '-', '0');
pscl_wD4 = strrep(pscl_wD4, 'V', '4');
% 转换为数值数组
pscl_wD4 = str2double(pscl_wD4);
NEWSUID_wD4 = NEWSUID(2:end,1);

[wD5,NEWSUID,~] =xlsread('F:\Google下载\wise_30sec_v1\wise_30sec_v1\导出实验\HW30s_wD5_to_excel.xls','HW30s_wD5_to_excel');
NEWSUID = string(NEWSUID);
pscl_wD5 = NEWSUID(2:end,23);
% 多重替换
pscl_wD5 = strrep(pscl_wD5, 'M', '2');
pscl_wD5 = strrep(pscl_wD5, 'F', '3'); 
pscl_wD5 = strrep(pscl_wD5, 'C', '1');
pscl_wD5 = strrep(pscl_wD5, '-', '0');
pscl_wD5 = strrep(pscl_wD5, 'V', '4');
% 转换为数值数组
pscl_wD5 = str2double(pscl_wD5);
NEWSUID_wD5 = NEWSUID(2:end,1);

bd_value = zeros(size(latlon(:,1),1),50); % 行为点，列为深度增加
for i=1:size(latlon(:,1),1)
    NEWSUID = bdgis_value(i);
    % 如果为0，就是空值，跳出
    if NEWSUID ==0
        continue;
    end  
    % 搜索在第几行
    [row, ~] = find(NEWSUID == wise30(:,1));
    % 取出对应的code
    CODEname = CODE(row);
    % 生成NEWSUID索引
    ling = length(num2str(NEWSUID));
    switch ling
        case 1
            SUID = 'WD' + CODEname + '000000' + num2str(NEWSUID);
        case 2
            SUID = 'WD' + CODEname + '00000' + num2str(NEWSUID);
        case 3
            SUID = 'WD' + CODEname + '0000' + num2str(NEWSUID);
        case 4
            SUID = 'WD' + CODEname + '000' + num2str(NEWSUID);
        case 5
            SUID = 'WD' + CODEname + '00' + num2str(NEWSUID);
    end
    % 取出容重值
    [row, ~] = find(SUID == NEWSUID_wD1);
    % Sand SDTO 16
    bd_value(i,1) = wD1(row,16);
    bd_value(i,11) = wD2(row,16);
    bd_value(i,21) = wD3(row,16);
    bd_value(i,31) = wD4(row,16);
    bd_value(i,41) = wD5(row,16);
    % Silt STPC 18
    bd_value(i,2) = wD1(row,18);
    bd_value(i,12) = wD2(row,18);
    bd_value(i,22) = wD3(row,18);
    bd_value(i,32) = wD4(row,18);
    bd_value(i,42) = wD5(row,18);
    % Clay CLPC 20
    bd_value(i,3) = wD1(row,20);
    bd_value(i,13) = wD2(row,20);
    bd_value(i,23) = wD3(row,20);
    bd_value(i,33) = wD4(row,20);
    bd_value(i,43) = wD5(row,20);
    % texture PSCL
    bd_value(i,4) = pscl_wD1(row,1);
    bd_value(i,14) = pscl_wD2(row,1);
    bd_value(i,24) = pscl_wD3(row,1);
    bd_value(i,34) = pscl_wD4(row,1);
    bd_value(i,44) = pscl_wD5(row,1);
    % BD BULK 23
    bd_value(i,5) = wD1(row,23);
    bd_value(i,15) = wD2(row,23);
    bd_value(i,25) = wD3(row,23);
    bd_value(i,35) = wD4(row,23);
    bd_value(i,45) = wD5(row,23);
    % SOC_gkg ORGC 27
    bd_value(i,6) = wD1(row,27);
    bd_value(i,16) = wD2(row,27);
    bd_value(i,26) = wD3(row,27);
    bd_value(i,36) = wD4(row,27);
    bd_value(i,46) = wD5(row,27);
    % TN_gkg TOTN 29
    bd_value(i,7) = wD1(row,29);
    bd_value(i,17) = wD2(row,29);
    bd_value(i,27) = wD3(row,29);
    bd_value(i,37) = wD4(row,29);
    bd_value(i,47) = wD5(row,29);
    % CN CNrt 31
    bd_value(i,8) = wD1(row,31);
    bd_value(i,18) = wD2(row,31);
    bd_value(i,28) = wD3(row,31);
    bd_value(i,38) = wD4(row,31);
    bd_value(i,48) = wD5(row,31);
    % PH PHAQ 33
    bd_value(i,9) = wD1(row,33);
    bd_value(i,19) = wD2(row,33);
    bd_value(i,29) = wD3(row,33);
    bd_value(i,39) = wD4(row,33);
    bd_value(i,49) = wD5(row,33);
    % CEC ECEC 37
    bd_value(i,10) = wD1(row,37);
    bd_value(i,20) = wD2(row,37);
    bd_value(i,30) = wD3(row,37);
    bd_value(i,40) = wD4(row,37);
    bd_value(i,50) = wD5(row,37);
end

% 后来加上提取的I值
% 土壤容重数据,g cm-3，用于转换单位
BD1 = bd_value(:,5);
BD2 = bd_value(:,15);
BD3 = bd_value(:,25);
BD4 = bd_value(:,35);
BD5 = bd_value(:,45);
clear NEWSUID NEWSUID_wD1 NEWSUID_wD2 NEWSUID_wD3 NEWSUID_wD4 NEWSUID_wD5 pscl_wD1 pscl_wD2 pscl_wD3 pscl_wD4 pscl_wD5
clear wD1 wD2 wD3 wD4 wD5 SUID wise30 filename CODE CODEname bd_value bdgis_value
clear ling i

% 碳输入已知  g C m-2 y-1. 偏移系数10  
[NPP, R] = geotiffread('H:/数据库/NPP_NDVI_MODIS/土壤总NPP分配20cm间隔.tif');
[row, col] = geographicToDiscrete(R, latlon(:,1), latlon(:,2));  % 转换为栅格行列号
temp=NPP(:,:,1);
NPP1 = temp(sub2ind(size(temp), row, col)); %取出对应位置的值
temp=NPP(:,:,2);
NPP2 = temp(sub2ind(size(temp), row, col)); %取出对应位置的值
temp=NPP(:,:,3);
NPP3 = temp(sub2ind(size(temp), row, col)); %取出对应位置的值
temp=NPP(:,:,4);
NPP4 = temp(sub2ind(size(temp), row, col)); %取出对应位置的值
temp=NPP(:,:,5);
NPP5 = temp(sub2ind(size(temp), row, col)); %取出对应位置的值
clear NPP temp

NPP1 = double(NPP1);
NPP2 = double(NPP2);
NPP3 = double(NPP3);
NPP4 = double(NPP4);
NPP5 = double(NPP5);
NPP1(NPP1==32767) = NaN;
NPP2(NPP2==32767) = NaN;
NPP3(NPP3==32767) = NaN;
NPP4(NPP4==32767) = NaN;
NPP5(NPP5==32767) = NaN;

% g C m-2 y-1 --> g C kg-1 day-1
NPP1 = NPP1/10./(BD1*100*100*20/1000)/365;
NPP2 = NPP2/10./(BD2*100*100*20/1000)/365;
NPP3 = NPP3/10./(BD3*100*100*20/1000)/365;
NPP4 = NPP4/10./(BD4*100*100*20/1000)/365;
NPP5 = NPP5/10./(BD5*100*100*20/1000)/365;

I0_30 = mean(cat(3,NPP1,NPP1,NPP2), 3,'omitnan');
I30_100 = mean(cat(3,NPP2,NPP3,NPP4,NPP5,NPP3,NPP4,NPP5), 3,'omitnan');
clear NPP1 NPP2 NPP3 NPP4 NPP5

param2(:,12) = I0_30;
param2(:,13) = I30_100;
clear I0_30 I30_100

% 又补充了Qmax数据  0.5度 由原来的数据直接进行单位转换，别的每没变
[Qmax0_30, R] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\Qmax0-30cm.tif');
[Qmax30_100, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\Qmax30-100cm.tif');
[row, col] = geographicToDiscrete(R, param2(:,1), param2(:,2));  % 转换为栅格行列号
Qmax0_30 = Qmax0_30(sub2ind(size(Qmax0_30), row, col)); %取出对应位置的值
Qmax30_100 = Qmax30_100(sub2ind(size(Qmax30_100), row, col)); %取出对应位置的值
% Qmax = [Qmax0_30, Qmax30_100];% 两列 0-30 30-100
param2(:,14) = Qmax0_30;
param2(:,15) = Qmax30_100;
clear Qmax0_30 Qmax30_100

% POC 用原始数据 得将20cm间隔的做平均 
% 空值用-1填充, kg C m-2，0.5度
[POC20, R] = geotiffread('H:\BaiduSyncdisk\基金\青年基金\来自文献的残体汇编数据\POC_MAOC数据周正虎NC文章中\POC0-20cm.tif');
[POC40, ~] = geotiffread('H:\BaiduSyncdisk\基金\青年基金\来自文献的残体汇编数据\POC_MAOC数据周正虎NC文章中\POC20-40cm.tif');
[POC60, ~] = geotiffread('H:\BaiduSyncdisk\基金\青年基金\来自文献的残体汇编数据\POC_MAOC数据周正虎NC文章中\POC40-60cm.tif');
[POC80, ~] = geotiffread('H:\BaiduSyncdisk\基金\青年基金\来自文献的残体汇编数据\POC_MAOC数据周正虎NC文章中\POC60-80cm.tif');
[POC100, ~] = geotiffread('H:\BaiduSyncdisk\基金\青年基金\来自文献的残体汇编数据\POC_MAOC数据周正虎NC文章中\POC80-100cm.tif');
mask2 = (POC20<=0) | (POC40 <= 0) | (POC60 <= 0) | (POC80 <= 0) | (POC100 <= 0);
POC20(mask2) = NaN;
POC40(mask2) = NaN;
POC60(mask2) = NaN;
POC80(mask2) = NaN;
POC100(mask2) = NaN;
clear mask2
[row, col] = geographicToDiscrete(R, param2(:,1), param2(:,2));  % 转换为栅格行列号
POC20 = POC20(sub2ind(size(POC20), row, col)); %取出对应位置的值
POC40 = POC40(sub2ind(size(POC40), row, col)); %取出对应位置的值
POC60 = POC60(sub2ind(size(POC60), row, col)); %取出对应位置的值
POC80 = POC80(sub2ind(size(POC80), row, col)); %取出对应位置的值
POC100 = POC100(sub2ind(size(POC100), row, col)); %取出对应位置的值

% 使用容重进行单位转换
% kg C m2 --> g C kg-1
POC0_20 = POC20./(BD1*100*100*20/1000).*1000;
POC20_40 = POC40./(BD2*100*100*20/1000).*1000;
POC40_60 = POC60./(BD3*100*100*20/1000).*1000;
POC60_80 = POC80./(BD4*100*100*20/1000).*1000;
POC80_100 = POC100./(BD5*100*100*20/1000).*1000;
% 将20cm间隔的做平均 
POC0_30 = mean(cat(3,POC0_20,POC0_20,POC20_40), 3,'omitnan');
POC30_100 = mean(cat(3,POC20_40,POC40_60,POC60_80,POC80_100,POC40_60,POC60_80,POC80_100), 3,'omitnan');
clear POC0_20 POC20_40 POC40_60 POC60_80 POC80_100 POC20 POC40 POC60 POC80 POC100
param2(:,16) = POC0_30;
param2(:,17) = POC30_100;
clear POC0_30 POC30_100

% DOC  没有进行深度转换，中间深度小于30的用0-30的，中间深度大于30的用30-100的，然后用容重转换单位
% 0-30cm  30-100cm,   g C kg-1，0.5度 
[DOC, R] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\估计的DOC0-30cm_30-100cm.tif'); 
[row, col] = geographicToDiscrete(R, param2(:,1), param2(:,2));  % 转换为栅格行列号
DOC0_30=DOC(:,:,1);
DOC30_100=DOC(:,:,2);
DOC0_30 = DOC0_30(sub2ind(size(DOC0_30), row, col)); %取出对应位置的值
DOC30_100 = DOC30_100(sub2ind(size(DOC30_100), row, col)); %取出对应位置的值
param2(:,18) = DOC0_30;
param2(:,19) = DOC30_100;
clear DOC DOC0_30 DOC30_100


% MBC前面已计算
param2(:,20) = param2_mbc0_30;
param2(:,21) = param2_mbc30_100;
clear param2_mbc0_30 param2_mbc30_100


% MAOC  g C kg-1，0.5度
[MAOC0_30, R] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\MAOC0-30cm.tif');  
[MAOC30_100, ~] = geotiffread('H:\数据库\基于青年基金收集的数据\已知参数空间分布数据\MAOC30-100cm.tif');
[row, col] = geographicToDiscrete(R, param2(:,1), param2(:,2));  % 转换为栅格行列号
param2(:,22) = MAOC0_30(sub2ind(size(MAOC0_30), row, col)); %取出对应位置的值
param2(:,23) = MAOC30_100(sub2ind(size(MAOC30_100), row, col)); %取出对应位置的值
clear MAOC0_30 MAOC30_100

% BD
param2(:,24) = BD1;
param2(:,25) = BD2;
param2(:,26) = BD3;
param2(:,27) = BD4;
param2(:,28) = BD5;

% 1.lat 2.lon 3.kbd_0_30 4.kbd_30_100 5.kl_0_30 6.kl_30_100 7.kp_0_30 8.kp_30_100 9.kq_0_30 10.kq_30_100 
% 11.cue 12.I0_30 13.I30_100 14.Qmax0_30 15.Qmax30_100 16.POC0_30 17.POC30_100 18.DOC0_30 19.DOC30_100
% 20.mbc0_30 21.mbc30_100 22.MAOC0_30 23.MAOC30_100 24.BD1 25.BD2 26.BD3 27.BD4 28.BD5
Jia_site_param = param2(:,1:28);
clear BD1 BD2 BD3 BD4 BD5 col latlon param2 R row
save('H:\数据库\基于青年基金收集的数据\Manuscript\SOCTimeSeries-main\Jia_site_param.mat', 'Jia_site_param', '-v6');


%% 模型公式
% 参数：I kp kbd kq(klm kml) kl cue 
function out = Cmodel_vectorized(Cpool, para)
    % 输入说明：
    % Cpool - 120x5矩阵 [P,L,B,M,CO2pool]（每行一个时间步）
    % para - 120x9矩阵 [kp,kbd,kq,kl,cue,in,pi,kml,qmax]（每行对应时间步参数）
    
    % 提取各碳库（自动广播为120x1向量）
    P = Cpool(:,1);
    L = Cpool(:,2);
    B = Cpool(:,3);
    M = Cpool(:,4);
    
    % 参数解包（每列对应一个参数）
    pi_val = para(:,7);
    kp = para(:,1);
    kbd = para(:,2);
    kq = para(:,3);
    kml = para(:,8);
    klm = kml.*kq;  % 元素乘法
    qmax = para(:,9);
    kl = para(:,4);
    cue = para(:,5);
    in_val = para(:,6);
    
    % 向量化计算通量（.*表示元素乘法）
    Fpl = kp.*P;
    Fbm = kbd.*B;   % ,注意这里的kbd为kbd*mbc
    Fml = kml.*M./qmax;
    Flm = klm.*L.*(1-M./qmax);
    Flb = kl.*L;
    
    % 计算输出（120x4矩阵）
    out = [pi_val.*in_val - Fpl, ...
          (1-pi_val).*in_val + Fpl + Fbm + Fml - Flm - Flb, ...
          cue.*Flb - Fbm, ...
          Flm - Fml];
end


%% 稳态解析解公式
% 用公式求解稳态值
function out = wentai_vectorized(para)
    nan_mask = any(isnan(para), 2);       % 标记含NaN的行
    neg_mask = any(para< 0, 2);          % 标记含负值的行
    target_rows = nan_mask | neg_mask; % 同时满足条件的行逻辑索引
    para(target_rows, :) = NaN;
    % 输入说明：
    % para - 120x9矩阵 [kp,kbd,kq,kl,cue,in,pi,kml,qmax]（每行对应一组参数）
    
    % 参数解包（每列对应一个参数）
    pi_val = para(:,7);
    kp = para(:,1);
    kbd = para(:,2); % ,注意这里的kbd为kbd*mbc
    kq = para(:,3);
    kml = para(:,8);
    klm = kml.*kq;  % 元素乘法
    qmax = para(:,9);
    kl = para(:,4);
    cue = para(:,5);
    in_val = para(:,6);
    
    % 向量化计算（所有运算符改为元素级运算）
    % ,注意这里的kbd为kbd*mbc
%     B = sqrt((in_val./kbd).*cue./(1-cue));
    B = (in_val./kbd).*cue./(1-cue);
%     L = kbd.*B.^2./(cue.*kl);
    L = kbd.*B./(cue.*kl);
    P = pi_val.*in_val./kp;
    M = klm.*L.*qmax./(kml + klm.*L);
    
    % 输出120x4矩阵（每行对应一组参数的[P,L,B,M]）
    out = [P, L, B, M];
end

%% 计算去除异常值后的平均值
function [out2,out3,out4] = time_yichang_2d(data)
    % 二维数据异常值处理函数
    % 输入：data - 二维数组（如120行4列）
    % 输出：out - 2行4列数组，第一行为四分位数方法结果，第二行为中位数方法结果
    
    % 获取数据的维度
    [~, col] = size(data);
    
    % 初始化输出数组
    out2 = zeros(1, col);  % 中位数方法结果
    out3 = zeros(1, col);  % 中位数方法结果
    out4 = zeros(1, col);  % 中位数方法结果
	for ss = 1:col
        % 方法2：使用中位数方法检测并删除离群值
        clean_data = rmoutliers(squeeze(data(:,ss))', 'median');
        out2(1,ss) = mean(clean_data, 'omitnan');
        out3(1,ss) = std(clean_data, 'omitnan');
        
        n = size(clean_data(~isnan(clean_data)),2); % 计算样本大小
        t_alpha_2_n_minus_1 = tinv(1-0.05/2, n-1); % 查找t分布的临界值，对于95%置信水平使用0.975（因为tinv需要双侧的概率）
        out4(1,ss) = t_alpha_2_n_minus_1 * (out3(1,ss) / sqrt(n));
	end
        
end

%% 2. 定义核心计算函数（适配任意环境变量）
function [bin_centers, prob_source_weighted, prob_sink_weighted, var_min, var_max, bin_edges] = calc_env_prob(env_data, meanMAT, meanMAP, meanAI1, meanET0, meanNPP, POC, MAOC, TOC_HIST, TOC, n_bins)
    % 步骤1：生成统一有效掩码（过滤所有无效值）
    valid_mask = (~isnan(TOC)) & (~isnan(meanMAT)) & (meanMAT > -14) & (meanMAT < 29) & ...
        (~isnan(meanMAP)) & (meanMAP > 200) & (meanMAP < 3000) & (~isnan(meanAI1)) & ...
        (meanAI1 > 600/10000) & (meanAI1 < 29000/10000) & (meanET0 ~= -9) & (meanET0 > 160) & (meanET0 < 2800) & ...
        (meanNPP ~= 0) & (~isnan(meanNPP)) & (meanNPP > 0.006) & (meanNPP < 1.57) & ...
        (POC > 0) & (POC < 20) & ...
        (MAOC > 0) & (MAOC < 70) & ...
        (TOC_HIST > 0) & (TOC_HIST < 100);

    % 步骤2：展平为一维配对数据
    env_1D = env_data(valid_mask);
    SOC_1D = TOC(valid_mask);
    fprintf('有效配对样本数 = %d\n', length(env_1D));

    % 步骤3：分位数分箱
    var_min = min(env_1D);
    var_max = max(env_1D);
    bin_edges = quantile(env_1D, linspace(0, 1, n_bins+1));
    bin_edges = unique(bin_edges);
    n_bins = length(bin_edges) - 1;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

    % 步骤4：计算加权概率
    prob_source_weighted = zeros(1, n_bins);
    prob_sink_weighted = zeros(1, n_bins);
    for i = 1:n_bins
        bin_mask = (env_1D >= bin_edges(i)) & (env_1D < bin_edges(i+1));
        if i == n_bins
            bin_mask = (env_1D >= bin_edges(i)) & (env_1D <= bin_edges(i+1));
        end
        SOC_bin = SOC_1D(bin_mask);
        n_bin = sum(bin_mask);
        
        if n_bin == 0
            prob_source_weighted(i) = 0;
            prob_sink_weighted(i) = 0;
            continue;
        end

        % 加权概率核心公式
        sum_source_abs = sum(abs(SOC_bin(SOC_bin < 0)));
        sum_sink = sum(SOC_bin(SOC_bin > 0));
        sum_total_abs = sum(abs(SOC_bin));
        
        if sum_total_abs == 0
            prob_source_weighted(i) = 0;
            prob_sink_weighted(i) = 0;
        else
            prob_source_weighted(i) = sum_source_abs / sum_total_abs;
            prob_sink_weighted(i) = sum_sink / sum_total_abs;
        end
    end
end

%% 2. 定义核心计算函数（适配任意环境变量）
function [bin_centers, prob_source_weighted, prob_sink_weighted, var_min, var_max, bin_edges] = calc_env_prob2(env_data, meanMAT, meanMAP, meanAI1, meanET0, meanNPP, TOC, n_bins)
    % 步骤1：生成统一有效掩码（过滤所有无效值）
    valid_mask = (~isnan(TOC)) & (~isnan(meanMAT)) & (meanMAT > -14) & (meanMAT < 29) & ...
        (~isnan(meanMAP)) & (meanMAP > 200) & (meanMAP < 3000) & (~isnan(meanAI1)) & ...
        (meanAI1 > 600/10000) & (meanAI1 < 29000/10000) & (meanET0 ~= -9) & (meanET0 > 160) & (meanET0 < 2800) & ...
        (meanNPP ~= 0) & (~isnan(meanNPP)) & (meanNPP > 0.006) & (meanNPP < 1.57);

    % 步骤2：展平为一维配对数据
    env_1D = env_data(valid_mask);
    SOC_1D = TOC(valid_mask);
    fprintf('有效配对样本数 = %d\n', length(env_1D));

    % 步骤3：分位数分箱
    var_min = min(env_1D);
    var_max = max(env_1D);
    bin_edges = quantile(env_1D, linspace(0, 1, n_bins+1));
    bin_edges = unique(bin_edges);
    n_bins = length(bin_edges) - 1;
    bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;

    % 步骤4：计算加权概率
    prob_source_weighted = zeros(1, n_bins);
    prob_sink_weighted = zeros(1, n_bins);
    for i = 1:n_bins
        bin_mask = (env_1D >= bin_edges(i)) & (env_1D < bin_edges(i+1));
        if i == n_bins
            bin_mask = (env_1D >= bin_edges(i)) & (env_1D <= bin_edges(i+1));
        end
        SOC_bin = SOC_1D(bin_mask);
        n_bin = sum(bin_mask);
        
        if n_bin == 0
            prob_source_weighted(i) = 0;
            prob_sink_weighted(i) = 0;
            continue;
        end

        % 加权概率核心公式
        sum_source_abs = sum(abs(SOC_bin(SOC_bin < 0)));
        sum_sink = sum(SOC_bin(SOC_bin > 0));
        sum_total_abs = sum(abs(SOC_bin));
        
        if sum_total_abs == 0
            prob_source_weighted(i) = 0;
            prob_sink_weighted(i) = 0;
        else
            prob_source_weighted(i) = sum_source_abs / sum_total_abs;
            prob_sink_weighted(i) = sum_sink / sum_total_abs;
        end
    end
end

function vec_out = iqr_remove_outlier(vec_in)
    q1 = prctile(vec_in,25);
    q3 = prctile(vec_in,75);
    iqr_val = q3 - q1;
    low = q1 - 1.5*iqr_val;
    high = q3 +1.5*iqr_val;
    vec_out = vec_in(vec_in >= low & vec_in <= high);
end