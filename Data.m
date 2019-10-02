fid = fopen('Group2_Case1.txt', 'rt');
data_cell = textscan(fid, '%f', 'CommentStyle', {'File', 'Data  Position     Value'}, 'CollectData', 1);
fclose(fid);
data = data_cell{1};