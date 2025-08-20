function opts = set_default(opts, name, val)
if ~isfield(opts,name) || isempty(opts.(name))
    opts.(name) = val;
end
end
