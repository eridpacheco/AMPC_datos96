function dt = datos96_parse_date_input(inputDate, refTime)
%DATOS96_PARSE_DATE_INPUT Convert date input to datetime with timezone.

if nargin < 2 || isempty(refTime)
    refTime = datetime.empty(0, 1);
end

if isempty(inputDate)
    dt = [];
    return;
end

if isdatetime(inputDate)
    dt = inputDate;
elseif isnumeric(inputDate)
    dt = datetime(inputDate, 'ConvertFrom', 'datenum');
elseif isstring(inputDate) || ischar(inputDate)
    dt = datetime(inputDate);
else
    error('datos96_parse_date_input:InvalidDateInput', ...
        'Formato de fecha no soportado: %s.', class(inputDate));
end

if isdatetime(refTime) && ~isempty(refTime)
    refTz = refTime.TimeZone;
    if isempty(dt.TimeZone) && ~isempty(refTz)
        dt.TimeZone = refTz;
    end
end
end