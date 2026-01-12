files = dir(fullfile('ggtrans*.mat'));
allData = [];
for i = 1:length(files)
    filePath = fullfile(files(i).name);
    data = load(filePath, 'resultss'); 
    allData = [allData; data.resultss(:)]; 
end
uniqueData = unique(allData);
save(fullfile( 'combinedData.mat'), 'uniqueData'); 
