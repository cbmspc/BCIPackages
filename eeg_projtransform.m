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



