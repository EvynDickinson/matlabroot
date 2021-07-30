function f(handle,fieldname)

% this function deletes the field given by fieldname from the 'UserData'
% field of handle

% make sure h is a handle
if ~ishandle(handle)
  error('Invalid handle');
end

% make sure h has a UserData field
props = get(handle);
if ~isfield(props,'UserData')
  error('Object has no UserData field'); 
end

% delete the field
ud = get(handle,'UserData');
rmfield(ud,fieldname);
