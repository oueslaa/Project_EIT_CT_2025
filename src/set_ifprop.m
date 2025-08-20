function set_ifprop(obj, prop, val)
try
    if isprop(obj,prop)
        obj.(prop) = val;
    end
catch
end
end
