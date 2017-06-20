function frame = framesFromMSER(mser)
    location = single(round(mser.Location));
    sz = mser.Axes;
    sz = max(1.6,1/4*sqrt(sz(:,1).*sz(:,2)));
    sz = single(sz);
    orientation = single(mser.Orientation);
    %orientation = zeros(length(location),1);
    frame=[location, sz, orientation];
end
