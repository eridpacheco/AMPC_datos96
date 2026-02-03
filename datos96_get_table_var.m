function v = datos96_get_table_var(tbl, name)
%DATOS96_GET_TABLE_VAR Return table variable or NaNs if missing.

if ismember(name, tbl.Properties.VariableNames)
    v = tbl.(name);
    if iscell(v)
        v = cellfun(@(x) double(x), v);
    end
    if islogical(v)
        v = double(v);
    end
    if iscategorical(v)
        v = double(v);
    end
else
    v = nan(height(tbl), 1);
end
end