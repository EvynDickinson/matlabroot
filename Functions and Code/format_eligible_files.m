
function FileNames = format_eligible_files(eligible_files,filler_string)
% FileNames = format_eligible_files(eligible_files,filler_string)
% filler_string default:  ' -- '


if nargin <2
    filler_string = ' -- ';
end

LOC = false(size(eligible_files));
for i = 1:size(eligible_files,2)
try LOC(:,i) = cellfun(@isnan,eligible_files(:,i));
catch
end
end
c = cellfun(@string,eligible_files);
c(LOC) = filler_string;
FileNames = join(c);