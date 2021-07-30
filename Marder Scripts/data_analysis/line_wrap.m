function line_hs = f(x,y,ylims,varargin)

% convert x and y to row vectors if they're not already
if ndims(x)==2 && size(x,2)==1
  x=x';
end
if ndims(y)==2 && size(y,2)==1
  y=y';
end

y_lo=ylims(1); y_hi=ylims(2);
y_span=y_hi-y_lo;
y_shift=y-y_lo;
y_shift_phase=y_shift/y_span;
y_shift_n=floor(y_shift_phase);
y_shift_normed_phase=y_shift_phase-y_shift_n;
i_breaks=find(diff(y_shift_n)~=0);
n_breaks=length(i_breaks);
if n_breaks==0
  line_hs=line(x,y,varargin{:});
else
  n_lines=n_breaks+1;
  y_shift_n_line=[y_shift_n(i_breaks) y_shift_n(end)];
  x_break=zeros(1,n_breaks);
  % calculate the x_breaks
  for i=1:n_breaks
    i_break=i_breaks(i);
    x_break(i)=interp1([y_shift_phase(i_break) y_shift_phase(i_break+1)],...
                       [x(i_break) x(i_break+1)],...
                       max(y_shift_n(i_break),y_shift_n(i_break+1)));
  end
  % make lines using the splices
  y_normed=y_span*y_shift_normed_phase+y_lo;
  y_normed_break_pre=y_span/2*(1+diff(y_shift_n_line))+y_lo;
  y_normed_break_post=y_span/2*(1-diff(y_shift_n_line))+y_lo; ;
  line_hs=zeros(1,n_lines);
  line_hs(1)=line([x(1:i_breaks(1)) x_break(1)],...
                  [y_normed(1:i_breaks(1)) y_normed_break_pre(1)],...
                  varargin{:});
  for i=2:n_breaks
    % make a line ending in the i'th break
    line_hs(i)=...
      line([x_break(i-1) x(i_breaks(i-1)+1:i_breaks(i)) x_break(i)],...
           [y_normed_break_post(i-1) ...
            y_normed(i_breaks(i-1)+1:i_breaks(i)) ...
            y_normed_break_pre(i)],...
           varargin{:});
  end
  line_hs(end)=line([x_break(end) x(i_breaks(end)+1:end)],...
                    [y_normed_break_post(end) y_normed(i_breaks(end)+1:end)],...
                    varargin{:});
end
