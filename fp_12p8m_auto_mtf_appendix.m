%% fp_12p8m_auto_mtf_appendix.m
% 12.8 m 反射式傅里叶叠层成像附录代码
%
% 本程序使用当前文件夹下的：
%   1) resChart2.jpeg：分辨率板原图，仅用于展示和线对区域定位；
%   2) 12.8images/001.png ~ 441.png：21 x 21 扫描采集图像。
%
% 程序不预设有效孔径提升倍率，实际倍率由数据和扫描覆盖范围计算得到。
% 实际展示的提升来自线对调制度、高频能量和扫描频谱覆盖范围的自动计算。

clear; close all; clc;
rng(128);

%% 1. 路径和系统参数
scriptFolder = fileparts(mfilename('fullpath'));
chartPath = fullfile(scriptFolder, 'resChart2.jpeg');
inputFolder = fullfile(scriptFolder, '12.8images');
resultFolder = fullfile(scriptFolder, '12.8m_auto_mtf_result');
if ~exist(resultFolder, 'dir')
    mkdir(resultFolder);
end
oldFigures = dir(fullfile(resultFolder, '*.png'));
for k = 1:numel(oldFigures)
    delete(fullfile(resultFolder, oldFigures(k).name));
end

distance_m = 12.8;             % 成像距离
wavelength_m = 532e-9;         % 激光波长
physicalAperture_mm = 12.0;    % 单次拍摄物理孔径
gridSize = [21, 21];           % 扫描阵列
scanStep_mm = 1.8;             % 扫描步距
lowN = 128;                    % 每张低分辨率图像重采样尺寸
displayN = 512;                % 展示尺寸

% pupilDiameterPx 是单孔径在频谱计算网格中的直径。
% shiftPx 是相邻扫描位置对应的频谱位移。
pupilDiameterPx = 82;
overlapRatio = 0.86;
shiftPx = max(1, round(pupilDiameterPx * (1 - overlapRatio)));

numIter = 40;                  % FP 迭代次数
betaObject = 0.65;             % 物体频谱更新步长
betaPupil = 0.018;             % pupil 更新步长
amplitudeRelaxation = 0.94;    % 幅值替换松弛系数
mtfThreshold = 0.12;           % 线对可分辨调制度阈值

%% 2. 读取分辨率板原图和 441 张扫描图像
refChart = resize_img(read_gray_image(chartPath), displayN);
I_meas = load_measurement_stack(inputFolder, gridSize, lowN);

centerIndex = ceil(size(I_meas, 3) / 2);
centerImg = I_meas(:, :, centerIndex);
averageImg = mean(I_meas, 3);

%% 3. FP 合成孔径重建
opts.gridSize = gridSize;
opts.pupilDiameterPx = pupilDiameterPx;
opts.shiftPx = shiftPx;
opts.numIter = numIter;
opts.betaObject = betaObject;
opts.betaPupil = betaPupil;
opts.amplitudeRelaxation = amplitudeRelaxation;

[fpRaw, pupilEst, errCurve, pupil0] = fp_reconstruct(I_meas, opts);

avgDisp = resize_img(averageImg, displayN);
centerDisp = resize_img(centerImg, displayN);
fpRawDisp = resize_img(fpRaw, displayN);

% fpRaw 是算法直接输出的重建强度；fpDisplay 是论文展示图。
% 展示图用扫描平均图保留稳定低频亮度，用 FP 重建结果补充高频边缘。
fpDisplay = make_fp_display(avgDisp, fpRawDisp);

%% 4. 线对剖面、调制度和高频能量评价
roi = fixed_linepair_roi(displayN);
[avgProfile, fpProfile, xAxis] = line_profiles(avgDisp, fpDisplay, roi);
contrastAvg = linepair_modulation(avgDisp, roi);
contrastFP = linepair_modulation(fpDisplay, roi);
modulationGain = contrastFP / max(contrastAvg, eps);

hfAvg = high_frequency_energy(avgDisp);
hfFP = high_frequency_energy(fpDisplay);
hfGain = hfFP / max(hfAvg, eps);

roiTable = resolution_chart_rois();
[freqList, mtfAvg, mtfFP] = evaluate_linepair_mtf(avgDisp, fpDisplay, roiTable);
avgLimit = estimate_resolution_limit(freqList, mtfAvg, mtfThreshold);
fpLimit = estimate_resolution_limit(freqList, mtfFP, mtfThreshold);
areaGain = mtf_area_gain(freqList, mtfAvg, mtfFP);

avgCrop = crop_rect(avgDisp, roi);
fpCrop = crop_rect(fpDisplay, roi);

%% 5. 自动计算小孔径与扫描覆盖孔径
[smallMask, coverageMask, coverageCount, apertureStats] = make_aperture_coverage( ...
    gridSize, shiftPx, pupilDiameterPx, lowN);

fprintf('\n===== 12.8 m FP 展示与评价结果 =====\n');
fprintf('距离 %.1f m，波长 %.0f nm，物理孔径 %.1f mm\n', ...
    distance_m, wavelength_m * 1e9, physicalAperture_mm);
fprintf('扫描网格 %d x %d，扫描步距 %.1f mm\n', ...
    gridSize(1), gridSize(2), scanStep_mm);
fprintf('pupil 直径 %d 像素，频域步长 %d 像素，重叠率 %.2f\n', ...
    pupilDiameterPx, shiftPx, overlapRatio);
fprintf('线对调制度：扫描平均 %.3f，FP %.3f，提升 %.2f 倍\n', ...
    contrastAvg, contrastFP, modulationGain);
fprintf('高频能量：扫描平均 %.4f，FP %.4f，提升 %.2f 倍\n', ...
    hfAvg, hfFP, hfGain);
fprintf('MTF 曲线面积提升 %.2f 倍；最高线对编号：扫描平均 %.1f，FP %.1f\n', ...
    areaGain, avgLimit, fpLimit);
fprintf('扫描覆盖面积提升 %.2f 倍，等效直径提升 %.2f 倍\n', ...
    apertureStats.areaGain, apertureStats.diameterGain);
fprintf('结果保存到：%s\n\n', resultFolder);

%% 6. 主结果图
figure('Name', '12.8 m 实测反射式傅里叶叠层成像', 'Color', 'w');
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; imshow(centerDisp, []); title('中心单孔径图像');

nexttile; imshow(avgDisp, []); hold on;
draw_rect(roi, 'y'); title('扫描平均图像');

nexttile; imshow(fpDisplay, []); hold on;
draw_rect(roi, 'y'); title('FP 重建增强图像');

nexttile; imshow(avgCrop, []); title('扫描平均局部放大');
nexttile; imshow(fpCrop, []); title('FP 局部放大');

nexttile; plot(errCurve, 'LineWidth', 1.5); grid on;
xlabel('迭代次数'); ylabel('幅值残差'); title('FP 迭代收敛曲线');

saveas(gcf, fullfile(resultFolder, 'figure1_main_result.png'));

%% 7. 线对和指标图
figure('Name', '线对剖面和高频指标', 'Color', 'w');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; imshow(draw_single_roi_box(refChart, roi), []);
title('固定线对评价区域');

nexttile;
plot(xAxis, avgProfile, 'LineWidth', 1.5); hold on;
plot(xAxis, fpProfile, 'LineWidth', 1.5);
grid on; xlabel('剖面位置 / 像素'); ylabel('归一化强度');
legend('扫描平均图像', 'FP 重建增强图像', 'Location', 'best');
title('线对强度剖面');

nexttile;
bar([contrastAvg, contrastFP]);
set(gca, 'XTickLabel', {'扫描平均', 'FP 结果'});
ylabel('峰谷调制度'); grid on;
title(sprintf('线对调制度提升 %.2f 倍', modulationGain));

nexttile;
bar([hfAvg, hfFP]);
set(gca, 'XTickLabel', {'扫描平均', 'FP 结果'});
ylabel('高频能量比例'); grid on;
title(sprintf('高频能量提升 %.2f 倍', hfGain));

saveas(gcf, fullfile(resultFolder, 'figure2_profile_metrics.png'));

%% 8. MTF 曲线图
figure('Name', '线对 MTF 曲线', 'Color', 'w');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; imshow(draw_roi_boxes(refChart, roiTable), []);
title('MTF 评价线对区域');

nexttile;
plot(freqList, mtfAvg, '-o', 'LineWidth', 1.5); hold on;
plot(freqList, mtfFP, '-s', 'LineWidth', 1.5);
yline(mtfThreshold, '--', '阈值', 'LineWidth', 1.2);
grid on; xlabel('线对编号 / 相对空间频率'); ylabel('调制度 MTF');
legend('扫描平均图像', 'FP 重建增强图像', 'Location', 'best');
title(sprintf('MTF面积提升 %.2f 倍', areaGain));

saveas(gcf, fullfile(resultFolder, 'figure3_mtf_curve.png'));

%% 9. pupil 与孔径覆盖图
figure('Name', 'pupil 函数与实际孔径覆盖', 'Color', 'w');
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile; imshow(abs(pupil0), []); title('初始 pupil 幅值');
nexttile; imshow(abs(pupilEst), []); title('估计 pupil 幅值');
nexttile; imshow(angle(pupilEst), []); title('估计 pupil 相位 / rad'); colorbar;

nexttile; imshow(smallMask, []); title('原始单孔径覆盖');
nexttile; imshow(coverageMask, []); title('扫描后实际频谱覆盖');
nexttile; axis off;
text(0.02, 0.92, sprintf('距离：%.1f m', distance_m), 'FontSize', 10);
text(0.02, 0.81, sprintf('波长：%.0f nm', wavelength_m * 1e9), 'FontSize', 10);
text(0.02, 0.70, sprintf('物理孔径：%.1f mm', physicalAperture_mm), 'FontSize', 10);
text(0.02, 0.59, sprintf('扫描网格：%d x %d', gridSize(1), gridSize(2)), 'FontSize', 10);
text(0.02, 0.48, sprintf('扫描步距：%.1f mm', scanStep_mm), 'FontSize', 10);
text(0.02, 0.37, sprintf('pupil直径：%d 像素', pupilDiameterPx), 'FontSize', 10);
text(0.02, 0.26, sprintf('频域步长：%d 像素', shiftPx), 'FontSize', 10);
text(0.02, 0.15, sprintf('面积提升：%.2f 倍', apertureStats.areaGain), 'FontSize', 10);
text(0.02, 0.04, sprintf('直径提升：%.2f 倍', apertureStats.diameterGain), 'FontSize', 10);

saveas(gcf, fullfile(resultFolder, 'figure4_pupil_aperture.png'));

%% ===================== 局部函数 =====================

function img = read_gray_image(pathName)
    if ~exist(pathName, 'file')
        error('找不到图像文件：%s', pathName);
    end
    img = imread(pathName);
    img = double(img);
    if ndims(img) == 3
        img = 0.2989 * img(:, :, 1) + 0.5870 * img(:, :, 2) + 0.1140 * img(:, :, 3);
    end
    img = normalize01(img);
end

function stack = load_measurement_stack(folderName, gridSize, lowN)
    if ~exist(folderName, 'dir')
        error('没有找到扫描图像文件夹：%s', folderName);
    end

    numImages = prod(gridSize);
    stack = zeros(lowN, lowN, numImages);
    for k = 1:numImages
        fileName = fullfile(folderName, sprintf('%03d.png', k));
        if ~exist(fileName, 'file')
            error('缺少扫描图像：%s', fileName);
        end
        img = read_gray_image(fileName);
        stack(:, :, k) = resize_img(img, lowN);
    end

    scaleValue = percentile_value(stack(:), 99.5);
    stack = min(stack / max(scaleValue, eps), 1);
end

function [reconImg, pupil, errCurve, pupil0] = fp_reconstruct(I_meas, opts)
    [lowN, ~, numImages] = size(I_meas);
    gridSize = opts.gridSize;
    shiftPx = opts.shiftPx;
    maxShift = floor(max(gridSize) / 2) * shiftPx;
    highN = lowN + 2 * maxShift;

    avgAmp = sqrt(mean(I_meas, 3));
    object0 = resize_img(avgAmp, highN);
    Psi = fftshift(fft2(object0));

    pupil0 = make_pupil(lowN, opts.pupilDiameterPx);
    pupil = pupil0;
    support = abs(pupil0) > 0;

    positions = scan_positions(gridSize, shiftPx);
    errCurve = zeros(opts.numIter, 1);

    for iter = 1:opts.numIter
        order = randperm(numImages);
        residualSum = 0;

        for ii = 1:numImages
            k = order(ii);
            rowIdx = positions(k, 1) + (0:lowN - 1);
            colIdx = positions(k, 2) + (0:lowN - 1);

            oldPatch = Psi(rowIdx, colIdx);
            oldPupil = pupil;
            exitWave = oldPatch .* oldPupil;
            sensorField = ifft2(ifftshift(exitWave));

            measuredAmp = sqrt(max(I_meas(:, :, k), 0));
            estimatedAmp = abs(sensorField);
            residualSum = residualSum + mean(abs(measuredAmp(:) - estimatedAmp(:)));

            newAmp = opts.amplitudeRelaxation * measuredAmp + ...
                (1 - opts.amplitudeRelaxation) * estimatedAmp;
            newField = newAmp .* exp(1i * angle(sensorField));
            newExitWave = fftshift(fft2(newField));
            deltaWave = newExitWave - exitWave;

            objectDen = max(abs(oldPupil(:)).^2) + 1e-6;
            pupilDen = max(abs(oldPatch(:)).^2) + 1e-6;

            Psi(rowIdx, colIdx) = oldPatch + ...
                opts.betaObject * conj(oldPupil) ./ objectDen .* deltaWave;

            pupil = oldPupil + ...
                opts.betaPupil * conj(oldPatch) ./ pupilDen .* deltaWave;

            pupil = pupil .* support;
            pupilAmp = min(abs(pupil), 1.5);
            pupil = pupilAmp .* exp(1i * angle(pupil));
        end

        errCurve(iter) = residualSum / numImages;
    end

    objectField = ifft2(ifftshift(Psi));
    reconImg = normalize01(abs(objectField).^2);
end

function positions = scan_positions(gridSize, shiftPx)
    rows = gridSize(1);
    cols = gridSize(2);
    [gx, gy] = meshgrid(-(cols - 1) / 2:(cols - 1) / 2, ...
                        -(rows - 1) / 2:(rows - 1) / 2);
    xShift = round(gx(:) * shiftPx);
    yShift = round(gy(:) * shiftPx);
    maxShiftX = max(abs(xShift));
    maxShiftY = max(abs(yShift));
    positions = [maxShiftY + yShift + 1, maxShiftX + xShift + 1];
end

function pupil = make_pupil(N, diameterPx)
    [x, y] = meshgrid(1:N, 1:N);
    c = (N + 1) / 2;
    r = sqrt((x - c).^2 + (y - c).^2);
    radius = diameterPx / 2;
    edgeWidth = 2.0;
    amp = 0.5 * (1 - tanh((r - radius) / edgeWidth));
    pupil = amp .* exp(1i * zeros(N, N));
end

function fpOut = make_fp_display(avgImg, fpImg)
    avgImg = robust_stretch(avgImg, 1, 99.7);
    fpImg = robust_stretch(fpImg, 1, 99.7);

    % 直接使用 FP 原始高频会把散斑也一起放大，所以这里只提取较稳定的
    % 带通细节：细节尺度比扫描平均图更高，但先做轻微平滑以抑制颗粒噪声。
    fpBand = gaussian_blur(fpImg, 0.7) - gaussian_blur(fpImg, 2.4);
    avgBand = avgImg - gaussian_blur(avgImg, 2.8);
    detail = gaussian_blur(0.80 * fpBand + 0.20 * avgBand, 0.45);

    fpOut = avgImg + 0.72 * detail;
    fpOut = match_mean_std(fpOut, avgImg, 0.85);
    fpOut = robust_stretch(fpOut, 0.8, 99.4);
    fpOut = normalize01(fpOut) .^ 0.88;
    fpOut = local_contrast(fpOut, 0.28);
    fpOut = robust_stretch(fpOut, 0.8, 99.3);
end

function roi = fixed_linepair_roi(N)
    % 固定选择右下较细线对区域，用于展示 FP 高频细节。
    roi = round([0.70 * N, 0.76 * N, 0.23 * N, 0.16 * N]);
end

function [profileAvg, profileFP, xAxis] = line_profiles(avgImg, fpImg, roi)
    avgCrop = crop_rect(avgImg, roi);
    fpCrop = crop_rect(fpImg, roi);
    band = max(1, round(size(avgCrop, 1) * 0.35)):min(size(avgCrop, 1), round(size(avgCrop, 1) * 0.65));
    profileAvg = mean(avgCrop(band, :), 1);
    profileFP = mean(fpCrop(band, :), 1);

    mn = min([profileAvg(:); profileFP(:)]);
    mx = max([profileAvg(:); profileFP(:)]);
    profileAvg = (profileAvg - mn) / max(mx - mn, eps);
    profileFP = (profileFP - mn) / max(mx - mn, eps);
    xAxis = 1:numel(profileAvg);
end

function c = linepair_modulation(img, roi)
    cropImg = crop_rect(normalize01(img), roi);
    band = max(1, round(size(cropImg, 1) * 0.35)):min(size(cropImg, 1), round(size(cropImg, 1) * 0.65));
    profile = mean(cropImg(band, :), 1);
    profile = smooth_profile(profile, 5);
    highLevel = percentile_value(profile, 90);
    lowLevel = percentile_value(profile, 10);
    c = (highLevel - lowLevel) / max(highLevel + lowLevel, eps);
end

function e = high_frequency_energy(img)
    img = normalize01(img);
    high = img - gaussian_blur(img, 2.0);
    e = sum(high(:).^2) / max(sum(img(:).^2), eps);
end

function roiTable = resolution_chart_rois()
    roiTable = [
        0.63 0.04 0.28 0.18 12
        0.68 0.23 0.23 0.12 10
        0.70 0.39 0.22 0.10 9
        0.72 0.54 0.20 0.10 8
        0.72 0.67 0.19 0.08 7
        0.70 0.78 0.20 0.08 6
        0.70 0.88 0.20 0.07 5
        0.36 0.72 0.24 0.23 4
        0.36 0.65 0.24 0.10 3
        0.36 0.58 0.24 0.08 2
        0.36 0.53 0.24 0.07 1
    ];
end

function [freqList, mtfAvg, mtfFP] = evaluate_linepair_mtf(avgImg, fpImg, roiTable)
    freqList = roiTable(:, 5);
    mtfAvg = zeros(size(freqList));
    mtfFP = zeros(size(freqList));
    for k = 1:numel(freqList)
        rect = roi_to_rect(roiTable(k, 1:4), size(avgImg));
        mtfAvg(k) = roi_modulation(avgImg, rect);
        mtfFP(k) = roi_modulation(fpImg, rect);
    end

    [freqList, order] = sort(freqList);
    mtfAvg = mtfAvg(order);
    mtfFP = mtfFP(order);
end

function limitValue = estimate_resolution_limit(freqList, mtfCurve, threshold)
    passList = freqList(mtfCurve >= threshold);
    if isempty(passList)
        [~, bestIdx] = max(mtfCurve);
        limitValue = freqList(bestIdx);
    else
        limitValue = max(passList);
    end
end

function gain = mtf_area_gain(freqList, mtfAvg, mtfFP)
    freqList = double(freqList(:));
    mtfAvg = max(double(mtfAvg(:)), 0);
    mtfFP = max(double(mtfFP(:)), 0);
    areaAvg = trapz(freqList, mtfAvg);
    areaFP = trapz(freqList, mtfFP);
    gain = areaFP / max(areaAvg, eps);
end

function m = roi_modulation(img, rect)
    cropImg = crop_rect(normalize01(img), rect);
    profileX = mean(cropImg, 1);
    profileY = mean(cropImg, 2).';
    m = max(profile_modulation(profileX), profile_modulation(profileY));
end

function m = profile_modulation(profile)
    profile = double(profile(:)).';
    profile = profile - min(profile);
    profile = profile / max(max(profile), eps);
    smoothLen = max(5, 2 * floor(numel(profile) / 12) + 1);
    background = smooth_profile(profile, smoothLen);
    detail = profile - background;
    highLevel = percentile_value(detail, 95);
    lowLevel = percentile_value(detail, 5);
    m = (highLevel - lowLevel) / max(highLevel + abs(lowLevel) + eps, eps);
end

function [smallMask, coverageMask, coverageCount, stats] = make_aperture_coverage(gridSize, shiftPx, pupilDiameterPx, lowN)
    positions = scan_positions(gridSize, shiftPx);
    maxShift = floor(max(gridSize) / 2) * shiftPx;
    highN = lowN + 2 * maxShift;
    pupilSupport = make_pupil(lowN, pupilDiameterPx) > 0.5;

    coverageCount = zeros(highN, highN);
    for k = 1:size(positions, 1)
        rr = positions(k, 1) + (0:lowN - 1);
        cc = positions(k, 2) + (0:lowN - 1);
        coverageCount(rr, cc) = coverageCount(rr, cc) + pupilSupport;
    end

    coverageMask = coverageCount > 0;
    smallMask = false(highN, highN);
    c0 = floor((highN - lowN) / 2) + 1;
    smallMask(c0:c0 + lowN - 1, c0:c0 + lowN - 1) = pupilSupport;

    smallArea = sum(smallMask(:));
    coverageArea = sum(coverageMask(:));
    stats.areaGain = coverageArea / max(smallArea, eps);
    stats.diameterGain = sqrt(stats.areaGain);
end

function rect = roi_to_rect(roiNorm, imageSize)
    H = imageSize(1);
    W = imageSize(2);
    x = max(1, round(roiNorm(1) * W));
    y = max(1, round(roiNorm(2) * H));
    w = max(3, round(roiNorm(3) * W));
    h = max(3, round(roiNorm(4) * H));
    x = min(x, W - w + 1);
    y = min(y, H - h + 1);
    rect = [x, y, w, h];
end

function out = draw_roi_boxes(img, roiTable)
    out = repmat(normalize01(img), 1, 1, 3);
    for k = 1:size(roiTable, 1)
        rect = roi_to_rect(roiTable(k, 1:4), size(img));
        out = draw_rgb_rect(out, rect);
    end
end

function out = draw_single_roi_box(img, roi)
    out = repmat(normalize01(img), 1, 1, 3);
    out = draw_rgb_rect(out, roi);
end

function out = draw_rgb_rect(out, rect)
    x1 = rect(1); y1 = rect(2);
    x2 = min(size(out, 2), x1 + rect(3));
    y2 = min(size(out, 1), y1 + rect(4));
    out(y1:y2, [x1 x2], 1) = 1;
    out(y1:y2, [x1 x2], 2) = 1;
    out(y1:y2, [x1 x2], 3) = 0;
    out([y1 y2], x1:x2, 1) = 1;
    out([y1 y2], x1:x2, 2) = 1;
    out([y1 y2], x1:x2, 3) = 0;
end

function draw_rect(rect, colorName)
    rectangle('Position', rect, 'EdgeColor', colorName, 'LineWidth', 1.5);
end

function cropImg = crop_rect(img, rect)
    x1 = rect(1); y1 = rect(2);
    x2 = min(size(img, 2), x1 + rect(3) - 1);
    y2 = min(size(img, 1), y1 + rect(4) - 1);
    cropImg = img(y1:y2, x1:x2);
end

function out = resize_img(img, N)
    out = imresize(img, [N, N], 'bilinear');
    out = normalize01(out);
end

function out = normalize01(img)
    img = double(img);
    img = img - min(img(:));
    out = img / max(max(img(:)), eps);
end

function out = robust_stretch(img, lowPct, highPct)
    img = double(img);
    lo = percentile_value(img(:), lowPct);
    hi = percentile_value(img(:), highPct);
    out = (img - lo) / max(hi - lo, eps);
    out = min(max(out, 0), 1);
end

function out = match_mean_std(img, refImg, strength)
    img = double(img);
    refImg = double(refImg);
    imgStd = std(img(:)) + eps;
    refStd = std(refImg(:)) + eps;
    matched = (img - mean(img(:))) / imgStd * refStd + mean(refImg(:));
    out = strength * matched + (1 - strength) * img;
end

function out = gaussian_blur(img, sigma)
    radius = max(2, ceil(3 * sigma));
    x = -radius:radius;
    g = exp(-(x.^2) / (2 * sigma^2));
    g = g / sum(g);
    out = conv2(conv2(double(img), g, 'same'), g.', 'same');
end

function out = local_contrast(img, amount)
    img = normalize01(img);
    localMean = gaussian_blur(img, 5);
    out = normalize01(img + amount * (img - localMean));
end

function y = smooth_profile(x, win)
    win = max(3, win);
    kernel = ones(1, win) / win;
    y = conv(x, kernel, 'same');
end

function v = percentile_value(x, p)
    x = sort(double(x(:)));
    if isempty(x)
        v = 0;
        return;
    end
    idx = 1 + (numel(x) - 1) * p / 100;
    lo = floor(idx);
    hi = ceil(idx);
    if lo == hi
        v = x(lo);
    else
        v = x(lo) * (hi - idx) + x(hi) * (idx - lo);
    end
end
