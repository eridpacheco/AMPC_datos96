function t = datos96_parse_timestamp(ts)
%DATOS96_PARSE_TIMESTAMP Robust conversion to datetime.

if isdatetime(ts)
    t = ts;
    return;
end

if isnumeric(ts)
    if all(ts > 1e12)
        t = datetime(ts / 1000, 'ConvertFrom', 'posixtime');
    elseif all(ts > 1e9)
        t = datetime(ts, 'ConvertFrom', 'posixtime');
    else
        t = datetime(ts, 'ConvertFrom', 'datenum');
    end
    return;
end

if isstring(ts) || iscellstr(ts) || ischar(ts)
    try
        t = datetime(ts);
    catch
        try
            t = datetime(ts, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        catch
            error('datos96_parse_timestamp:InvalidTimestamp', ...
                'No se pudo interpretar timestamp como datetime.');
        end
    end
    return;
end

error('datos96_parse_timestamp:UnsupportedType', ...
    'Tipo de timestamp no soportado: %s.', class(ts));
end