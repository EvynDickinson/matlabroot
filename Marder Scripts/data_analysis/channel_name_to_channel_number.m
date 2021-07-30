function channel_number=channel_name_to_channel_number(fid,channel_name)

% Determines a channel number given a channel name, for .smr files.

channel_list=SONChanList(fid);
channel_numbers=cell2mat({channel_list.number})';
channel_names={channel_list.title}';
channel_number=NaN;
for j=1:length(channel_names)
  if strcmp(channel_names{j},channel_name)
    channel_number=channel_numbers(j);
    break;
  end
end
if ~isfinite(channel_number)
  error('No channel by that name')
end
