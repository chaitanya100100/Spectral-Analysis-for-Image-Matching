function frame = framesFromMSER(mser)
% this function takes mser regions as input and returns sift frames so that
% they can be used as custom frame input in vl_sift function call.

    location = single(round(mser.Location));
    sz = mser.Axes;
    sz = max(1.6,1/4*sqrt(sz(:,1).*sz(:,2)));
    sz = single(sz);
    orientation = single(mser.Orientation);
    %orientation = zeros(length(location),1);
    frame=[location, sz, orientation];
end
