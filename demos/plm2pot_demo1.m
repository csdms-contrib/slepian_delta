%% PLM2POT DEMO1
function plm2pot_demo1
    % Make some surface mass in a region easy to plot
    [~, C, dels, dems, ~] = localization(18, 'samerica');
    C{1} = C{1} * 2.94; % Scale by 2.94 so the peak is roughly 50 kg/m^2 (like avg Greenland Melt)
    % Say the best function is our surface mass and convert it to geoid.
    lmcosiG = plm2pot([dels dems C{1}], [], [], [], 5);
    % Now put it back in surface density
    lmcosiSD = plm2pot(lmcosiG, [], [], [], 4);
    % Now make a figure of these 3 things
    ah = krijetem(subnum(3, 1));
    fig2print(gcf, 'tall')
    % Plot the spatial data in Mollweide projection
    axes(ah(1))
    [~, ~, ~] = plotplm([dels dems C{1}], [], [], 4, 1);
    kelicol
    t = title('Slepian fn from LOCALIZATION representing our surface density [kg/m^2]');
    %caxis([-13 1])
    cb = colorbar('hor');
    shrink(cb, 1.5, 2)
    movev(cb, -0.05)
    axes(cb)
    set(get(cb, 'xlabel'), 'string', 'Surface density [kg/m^2]');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(ah(2))
    % Plot in mm
    lmcosiG(:, 3:4) = lmcosiG(:, 3:4) * 1000;
    [~, ~, ~] = plotplm(lmcosiG, [], [], 4, 1);
    kelicol
    t = title('Panel 1 represented as Geoid [milimeters]');
    %caxis([-0.5 0.05])
    cb = colorbar('hor');
    shrink(cb, 1.5, 2)
    movev(cb, -0.05)
    axes(cb)
    set(get(cb, 'xlabel'), 'string', 'Geoid height variation [milimeters]');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    axes(ah(3))
    [~, ~, ~] = plotplm(lmcosiSD, [], [], 4, 1);
    kelicol
    t = title('Panel 2 converted back to surface density [kg/m^2]');
    %caxis([-13 1])
    cb = colorbar('hor');
    shrink(cb, 1.5, 2)
    movev(cb, -0.05)
    axes(cb)
    set(get(cb, 'xlabel'), 'string', 'Surface density [kg/m^2]');
end
