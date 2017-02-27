function CoordLookupTable = eeg_projtransform ()
% Creates a lookup table that transforms from a square grid to a grid
% projected from a hemisphere


% How many points along the radius (do not count the origin itself, but
% count the point landing on circumference) ?
Npts = 6;


k = 0;
for i = -Npts:Npts
    for j = -Npts:Npts
        k = k+1;
        CoordLookupTable.sqr(k,:) = [i,j];
    end
end

CoordLookupTable.proj = nan((Npts*2+1)^2,2);

for n = [-Npts+1:-1,1:Npts-1]
    % Center point
    r1 = n/Npts;
    t1 = pi/2;
    [x1,y1] = pol2cart(t1,r1);
    
    % Right point
    r2 = 1;
    t2 = pi/2/Npts*n;
    [x2,y2] = pol2cart(t2,r2);
    
    % Left point
    r3 = 1;
    t3 = pi/2/Npts*(2*Npts-n);
    [x3,y3] = pol2cart(t3,r3);
    
    % Find center of arching circle
    acircle = find_circle([x1,y1],[x2,y2],[x3,y3]);
    
    % Find intercepts of the two circles
    [xinter,yinter] = circcirc(0,0,1,acircle(1),acircle(2),acircle(3));

    % find which solution is right/left point
    ir = find(xinter>0);
    il = find(xinter<0);
    % re-define center to the center of arching circle
    acircle_intercept_right = [xinter(ir) yinter(ir)] - acircle(1:2);
    acircle_intercept_left = [xinter(il) yinter(il)] - acircle(1:2);
    % convert to polar
    [acircle_intercept_theta_right,acircle_intercept_radius] = cart2pol(acircle_intercept_right(1),acircle_intercept_right(2));
    [acircle_intercept_theta_left,acircle_intercept_radius] = cart2pol(acircle_intercept_left(1),acircle_intercept_left(2));

    % divide the arc into Npts*2+1 points (+1 to include the center point)
    arcdivision_theta = linspace(acircle_intercept_theta_left,acircle_intercept_theta_right,Npts*2+1);
    
    % Now walk on the arc of the arching circle that are inside the origin
    % circle. Always walk from left to right of the origin circle
    for i = 1:Npts*2+1
        % convert to cartesian, re-center to origin circle
        [x,y] = pol2cart(arcdivision_theta(i),acircle_intercept_radius);
        x = x+acircle(1);
        y = y+acircle(2);
        % This is the origin coordinate of this particular point
        CoordLookupTable.proj(eeg_findcoord(CoordLookupTable.sqr,[i-1-Npts,n]),:) = [x,y];
    end    
end

% When n=0, there is no arching circle
for j = -Npts:Npts
    CoordLookupTable.proj(eeg_findcoord(CoordLookupTable.sqr,[j,0]),:) = [j/Npts,0];
end
% When n=5, there is no arching circle
CoordLookupTable.proj(eeg_findcoord(CoordLookupTable.sqr,[0,Npts]),:) = [0,1];
CoordLookupTable.proj(eeg_findcoord(CoordLookupTable.sqr,[0,-Npts]),:) = [0,-1];



% function [xout,yout]=circcirc(x1,y1,r1,x2,y2,r2)
% %CIRCCIRC  Intersections of circles in Cartesian plane
% %
% %  [xout,yout] = CIRCCIRC(x1,y1,r1,x2,y2,r2) finds the points
% %  of intersection (if any), given two circles, each defined by center
% %  and radius in x-y coordinates.  In general, two points are
% %  returned.  When the circles do not intersect or are identical,
% %  NaNs are returned.  When the two circles are tangent, two identical
% %  points are returned.  All inputs must be scalars.
% %
% %  See also LINECIRC.
% 
% % Copyright 1996-2007 The MathWorks, Inc.
% % Written by:  E. Brown, E. Byrns
% 
% assert(isscalar(x1) && isscalar(y1) && isscalar(r1) && ...
%        isscalar(x2) && isscalar(y2) && isscalar(r2), ...
%     ['map:' mfilename ':mapError'], 'Inputs must be scalars')
% 
% assert(isreal([x1,y1,r1,x2,y2,r2]), ...
%     ['map:' mfilename ':mapError'], 'inputs must be real')
% 
% assert(r1 > 0 && r2 > 0, ...
%     ['map:' mfilename ':mapError'], 'radius must be positive')
% 
% % Cartesian separation of the two circle centers
% 
% r3=sqrt((x2-x1).^2+(y2-y1).^2);
% 
% indx1=find(r3>r1+r2);  % too far apart to intersect
% indx2=find(r2>r3+r1);  % circle one completely inside circle two
% indx3=find(r1>r3+r2);  % circle two completely inside circle one
% indx4=find((r3<10*eps)&(abs(r1-r2)<10*eps)); % circles identical
% indx=[indx1(:);indx2(:);indx3(:);indx4(:)];
% 
% anought=atan2((y2-y1),(x2-x1));
% 
% %Law of cosines
% 
% aone=acos(-((r2.^2-r1.^2-r3.^2)./(2*r1.*r3)));
% 
% alpha1=anought+aone;
% alpha2=anought-aone;
% 
% xout=[x1 x1]+[r1 r1].*cos([alpha1 alpha2]);
% yout=[y1 y1]+[r1 r1].*sin([alpha1 alpha2]);
% 
% % Replace complex results (no intersection or identical)
% % with NaNs.
% 
% if ~isempty(indx)
%     xout(indx,:) = NaN;    yout(indx,:) = NaN;
% end
