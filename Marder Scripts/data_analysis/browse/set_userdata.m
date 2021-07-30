function f(h,f,v)

%SET_USERDATA Set a field in the UserData of the specified object.
%   SET_USERDATA(H,F,V) Sets the field F of UserData from the 
%   object with handle H to V.

% set the appropriate field of the UserData structure to the 
% given value
ud = get(h,'UserData');
ud = setfield(ud,f,v);
set(h,'UserData',ud);
