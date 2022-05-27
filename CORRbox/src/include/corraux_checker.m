function output = corraux_checker(input)

cond1 = isstruct(input);
cond2 = isfield(input,'name');
cond3 = strcmp(input.name,'corr');
cond4 = (length(input.size)==3);
cond5 = all(input.size==size(input.data));

if (cond1&&cond2&&cond3&&cond4&&cond5)
    output = true;
else
    output = false;
end

end