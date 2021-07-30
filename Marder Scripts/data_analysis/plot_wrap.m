function line_hs = f(x,y,y_tol,varargin)

% convert x and y to row vectors if they're not already
if ndims(x)==2 && size(x,2)==1
  x=x';
end
if ndims(y)==2 && size(y,2)==1
  y=y';
end

i_break=find(abs(diff(y))>y_tol);
n_breaks=length(i_break);
x_break=zeros(1,n_breaks);
y_break_pre=zeros(1,n_breaks);
y_break_post=zeros(1,n_breaks);
% calculate the splices
x_pre=x(i_break);  x_post=x(i_break+1);
y_pre=y(i_break);  y_post=y(i_break+1);
for i=1:n_breaks
  if (y_post(i)>y_pre(i))
    y_post_new=y_post(i)-2*y_tol;
    x_break(i)=interp1([y_pre(i) y_post_new],[x_pre(i) x_post(i)],-y_tol);
    y_break_pre(i) =-y_tol;
    y_break_post(i)=+y_tol;    
  else
    y_post_new=y_post(i)+2*y_tol;
    x_break(i)=interp1([y_pre(i) y_post_new],[x_pre(i) x_post(i)],+y_tol);
    y_break_pre(i) =+y_tol;
    y_break_post(i)=-y_tol;    
  end
end
% make lines using the splices
line_hs=zeros(1,n_breaks+1);
line_hs(1)=line([x(1:i_break(1)) x_break(1)],...
                [y(1:i_break(1)) y_break_pre(1)],...
                varargin{:});
for i=2:n_breaks
  % make a line ending in the i'th break
  line_hs(i)=...
    line([x_break(i-1) x(i_break(i-1)+1:i_break(i)) x_break(i)],...
         [y_break_post(i-1) y(i_break(i-1)+1:i_break(i)) y_break_pre(i)],...
         varargin{:});
end
line_hs(end)=line([x_break(end) x(i_break(end)+1:end)],...
                  [y_break_post(end) y(i_break(end)+1:end)],...
                  varargin{:});

% make it like a plot
box on;