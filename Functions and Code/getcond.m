function [cond_label, labels] = getcond(cond)
%  [cond_label, labels] = getcond(cond)

load('conditionlabels');
labels = conditions;
cond_label = labels{cond};

end
