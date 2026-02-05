function datos96_add_irrigation_shading(figHandle, t, u_applied)
%DATOS96_ADD_IRRIGATION_SHADING Shade irrigation ON periods (control aplicado, uk_1).

if isempty(u_applied)
    return;
end
u_on = (u_applied > 0.5);
u_on(isnan(u_applied)) = false;
if all(~u_on)
    return;
end

figure(figHandle);
ax = gca;
yl = ylim(ax);

edges = diff([false; u_on(:); false]);
startIdx = find(edges == 1);
endIdx = find(edges == -1) - 1;

holdState = ishold(ax);
hold(ax, 'on');
for i = 1:numel(startIdx)
    x0 = t(startIdx(i));
    x1 = t(endIdx(i));
    patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], ...
        [0.9 0.9 0.9], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
if ~holdState
    hold(ax, 'off');
end

uistack(findobj(ax, 'Type', 'patch'), 'bottom');
end