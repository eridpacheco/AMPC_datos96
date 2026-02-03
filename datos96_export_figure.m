function datos96_export_figure(figHandle, pngPath, pdfPath)
%DATOS96_EXPORT_FIGURE Export figure to PNG (300 dpi) and PDF.

if ~exist(fileparts(pngPath), 'dir')
    mkdir(fileparts(pngPath));
end

try
    exportgraphics(figHandle, pngPath, 'Resolution', 300);
catch
    print(figHandle, pngPath, '-dpng', '-r300');
end

try
    exportgraphics(figHandle, pdfPath, 'ContentType', 'vector');
catch
    print(figHandle, pdfPath, '-dpdf');
end
end