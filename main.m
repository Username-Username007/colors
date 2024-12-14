%% Initialization
mkdir('gene')

%% Sensitivity of the Cones
close all;
red = makedist('normal', mu = 570, sigma = 50);
green = makedist('normal', mu = 540, sigma = 40);
blue = makedist('normal', mu = 450, sigma = 20);
    % the x should be wavelength

wavelen = linspace(380, 750);

fig = figure(color = 'w');
hold on;
plot(wavelen, color_signal('r', wavelen), 'r', linewidth = 2);
plot(wavelen, color_signal('g', wavelen), 'g', linewidth = 2);
plot(wavelen, color_signal('b', wavelen), 'b', linewidth = 2);
xlabel("Wavelength (nm)");
ylabel("Normalzied Signal");
exportgraphics(fig, 'gene/normalized_signal.png', resolution = 900, backgroundcolor = 'white');
%% Human Gamut
save_gif_path = "gene/gamut.gif";
close all;
wavelen = linspace(380, 750);
fig = figure(color = 'w');
hold on;
grid on;
rgb = [color_signal('r', wavelen)', color_signal('g', wavelen)', color_signal('b', wavelen)'];
% rgb = rgb ./ max(rgb, [], 1);
% scatter3(rgb(:, 1), rgb(:, 2), rgb(:, 3), 30, rgb ./ max(rgb,[], 2), 'filled', markeredgecolor = 'k');
for i = 1:(size(rgb, 1) - 1)
    plot3(rgb(i:i+1, 1), rgb(i:i+1, 2), rgb(i:i+1, 3), 'o-', color = rgb(i, :) / max(rgb(i, :), [], 'all'), LineWidth= 3);
end
xlabel('Long Cone (R)');
ylabel('Middle Cone (G)');
zlabel('Short Cone (B)');

for az = linspace(0, 360, 100)
    view(az, 20);
    fr = getframe(fig);
    [A, cmap] = rgb2ind(fr.cdata, 256);
    if az == 0
        imwrite(A, cmap, save_gif_path, 'gif', delaytime = .1, loopcount = inf);
    else
        imwrite(A, cmap, save_gif_path, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
    % pause(.1);
end

view([0, 0, 1])
exportgraphics(fig, 'gene/001.png', resolution = 900);
view([0, 1, 0])
exportgraphics(fig, 'gene/010.png', resolution = 900);
view([1, 0, 0])
exportgraphics(fig, 'gene/100.png', resolution = 900);
%% Spectrum
fig = figure(color = 'w');
hold on;
wavelen = linspace(350, 750, 500);
spectrum = [color_signal('r', wavelen)', color_signal('g', wavelen)', color_signal('b', wavelen)'];
spectrum = permute(repmat(spectrum, [1, 1, 50]), [3, 1, 2]);
% spectrum = spectrum ./ max(spectrum, [], 3);
imshow(spectrum, xData = [350, 750]);
axis on;
xlabel("Wavelength (nm)")
set(gca, 'YColor', 'none');
xlim([350, 750]);
exportgraphics(fig, 'gene/spectrum.png', resolution = 900);


%% CIE 1931 GIF
save_gif_path = 'gene/cie_1931.gif';
close all;
wavelen = linspace(380, 750);
fig = figure(color = 'w');
hold on;
rgb = [color_signal('r', wavelen)', color_signal('g', wavelen)', color_signal('b', wavelen)'];
for i = 1:(size(rgb, 1) - 1)
    plot3(rgb(i:i+1, 1), rgb(i:i+1, 2), rgb(i:i+1, 3), 'o-', color = rgb(i, :) / max(rgb(i, :), [], 'all'), LineWidth= 3, markersize = 3);
end

% Plot brightness line
for end_ponit_i = [25, 40, 50, 60]
    start_point = [0, 0, 0];
    end_point = rgb(end_ponit_i, :);
    mu = repmat(linspace(1, 0, 10)', [1, 3]);
    brightness_line = repmat(start_point, [length(mu), 1]) .* mu + (1-mu) .* repmat(end_point, [length(mu), 1]);
    scatter3(brightness_line(:, 1), brightness_line(:, 2), brightness_line(:, 3), 100, brightness_line ./ max(brightness_line(:)), 'filled');
end

view([1, .31, .6]);
xlabel('Long Cone (R)');
ylabel('Middle Cone (G)');
zlabel('Short Cone (B)');
grid on;

for az = linspace(0, 360, 200)
    view(az, 20);
    fr = getframe(fig);
    [A, cmap] = rgb2ind(fr.cdata, 256);
    if az == 0
        imwrite(A, cmap, save_gif_path, 'gif', delaytime = .1, loopcount = inf);
    else
        imwrite(A, cmap, save_gif_path, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
    % pause(.1);
end

%% CIE 1931 Plane
close all;
wavelen = linspace(380, 750);
fig = figure(color = 'w');
hold on;
grid on;
rgb = [color_signal('r', wavelen)', color_signal('g', wavelen)', color_signal('b', wavelen)'];
for i = 1:(size(rgb, 1) - 1)
    plot3(rgb(i:i+1, 1), rgb(i:i+1, 2), rgb(i:i+1, 3), 'o-', color = rgb(i, :) / max(rgb(i, :), [], 'all'), LineWidth= 3, markersize = 3);
end

% Plot the cutting surface
[X, Y] = meshgrid(linspace(0, .6));
Z = .6 - X - Y;
invalid_area = X<0 | Y<0 | Z<0;
X(invalid_area) = nan;
Y(invalid_area) = nan;
Z(invalid_area) = nan;
surf(X, Y, Z, edgecolor = 'none', facealpha = .9);

% Plot brightness line
for end_ponit_i = [25, 27, 30, 45]
    start_point = [0, 0, 0];
    end_point = rgb(end_ponit_i, :);
    mu = repmat(linspace(1, 0, 10)', [1, 3]);
    brightness_line = repmat(start_point, [length(mu), 1]) .* mu + (1-mu) .* repmat(end_point, [length(mu), 1]);
    scatter3(brightness_line(:, 1), brightness_line(:, 2), brightness_line(:, 3), 50, brightness_line ./ max(brightness_line(:)), 'filled');
    plot3(brightness_line(:, 1), brightness_line(:, 2), brightness_line(:, 3), color = brightness_line(end, :), linewidth = 2);
end

view([1, .1, .21]);
xlabel('Long Cone (R)');
ylabel('Middle Cone (G)');
zlabel('Short Cone (B)');
exportgraphics(fig, "gene/cie_1931_plane.png", resolution = 900);

%% Hyper green
save_gif_path = 'gene/hyper_green.gif';
close all;
wavelen = linspace(380, 750);
fig = figure(color = 'w');
hold on;
grid on;
rgb = [color_signal('r', wavelen)', color_signal('g', wavelen)', color_signal('b', wavelen)'];
% rgb = rgb ./ max(rgb, [], 1);
% scatter3(rgb(:, 1), rgb(:, 2), rgb(:, 3), 30, rgb ./ max(rgb,[], 2), 'filled', markeredgecolor = 'k');
for i = 1:(size(rgb, 1) - 1)
    plot3(rgb(i:i+1, 1), rgb(i:i+1, 2), rgb(i:i+1, 3), 'o-', color = rgb(i, :) / max(rgb(i, :), [], 'all'), LineWidth= 3);
end
xlabel('Long Cone (R)');
ylabel('Middle Cone (G)');
zlabel('Short Cone (B)');

scatter3(0, 1, 0, 400, 'g', 'filled');

grid on;

for az = linspace(0, 360, 200)
    view(az, 20);
    fr = getframe(fig);
    [A, cmap] = rgb2ind(fr.cdata, 256);
    if az == 0
        imwrite(A, cmap, save_gif_path, 'gif', delaytime = .1, loopcount = inf);
    else
        imwrite(A, cmap, save_gif_path, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
    % pause(.1);
end