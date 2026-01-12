
% 获取文件夹中所有的.mat文件
files = dir(fullfile('ggtrans*.mat'));

% 初始化一个空数组来存储所有数据
allData = [];

% 循环遍历每个文件
for i = 1:length(files)
    % 构建完整的文件路径
    filePath = fullfile(files(i).name);
    
  % 加载.mat文件中的数据
    data = load(filePath, 'resultss');  % 指定加载名为 resultss 的变量
    
    % 将加载的数据追加到allData数组中
    allData = [allData; data.resultss(:)];  % 确保数据是一维的
end

% 去除重复的数据
uniqueData = unique(allData);

% 保存合并后的数据到一个新的.mat文件中
save(fullfile( 'combinedData.mat'), 'uniqueData'); 