function [List] = getFileNumber(folder)
List = dir(folder);
List([1,2]) = [];
List = struct2table(List);
List(:,2:end) = [];

for n = 1:height(List)
    if height(List) > 1
        str1 = regexprep(char(List.name(n,1)),'[,;=]', ' ');
    else
        str1 = regexprep(char(List.name),'[,;=]', ' ');
    end
    str2 = regexprep(regexprep(str1,'[^- 0-9.eE(,)/]',''), ' \D* ',' ');
    str3 = regexprep(str2, {'\.\s','\E\s','\e\s','\s\E','\s\e'},' ');
    List.rpm(n,1) = str2double(str3);
end

if height(List) > 1
    [List.rpm, ListIdx] = sort(List.(2));
    List.name(1:end) = List.name(ListIdx);
end


