function roi = extract_roi(binary_img)
% ROIÃ·»°
stats = regionprops(binary_img, 'BoundingBox');
if ~isempty(stats)
    bbox = stats(1).BoundingBox;
    x = round(bbox(1));
    y = round(bbox(2));
    w = round(bbox(3));
    h = round(bbox(4));
    roi = binary_img(y:y+h-1, x:x+w-1);
else
    roi = binary_img;
end
end 