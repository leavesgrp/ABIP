function [] = write2csv(filename, title, results_file, format, type)
if nargin == 4
    type = "struct";
end
try
if type == "array"
    n_test = size(results_file,1);
    notitle_flag = false;
    if ~exist(filename,'file')
        notitle_flag = true;
    end
    fileID = fopen(filename, 'a');
    if fileID < 0
        disp('Fail to open file.');
    end

    % write title
    if notitle_flag
        for i = 1:length(title)-1
            fprintf(fileID,'%s,',title(i));
        end
        fprintf(fileID,'%s\n',title(end));
    end
    % write data
    for i = 1:n_test
        for j = 1:length(title)-1
            fprintf(fileID,strcat(format(j),','),results_file(i,j));
        end
        fprintf(fileID,strcat(format(length(title)),'\n'),results_file(i,length(title)));
    end

%     fclose(fileID);

elseif type == "struct"
    field = fieldnames(results_file);
    n_test = size(results_file,2);
    notitle_flag = false;
    if ~exist(filename,'file')
        notitle_flag = true;
    end
    fileID = fopen(filename, 'a');
    if fileID < 0
        disp('Fail to open file.');
    end

    % write title
    if notitle_flag
        for i = 1:length(title)-1
            fprintf(fileID,'%s,',title(i));
        end
        fprintf(fileID,'%s\n',title(end));
    end
    % write data
    for i = 1:n_test
        for j = 1:size(field,1)-1
            fprintf(fileID,strcat(format(j),','),results_file(i).(field{j}));
        end
        fprintf(fileID,strcat(format(length(title)),'\n'),results_file(i).(field{size(field,1)}));
    end
%     fclose(fileID);
end
catch
    fprintf("Fail to write file.\n");
end
fclose(fileID);
end