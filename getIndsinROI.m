function hits = getIndsinROI(subs,ROIs)
% for a given ROI, crop region 
numroi = size(ROIs,1);
hits = zeros(size(subs,1),1);
for iroi = 1:numroi
    roi = ROIs(iroi,:);
    rx = roi([1 4]);
    ry = roi(1+[1 4]);
    rz = roi(2+[1 4]);
    
    these = subs(:,1)>rx(1) & subs(:,1)<rx(2) &...
        subs(:,2)>ry(1) & subs(:,2)<ry(2)& ...
        subs(:,3)>rz(1) & subs(:,3)<rz(2);
    hits = hits | these;
end