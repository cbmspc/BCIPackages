% Find the center and radius of a circle, given three points that lie on
% the circle in cartesian coordinate system
% Method obtained from Yumnam Kirani Singh, 
% http://www.geocities.com/kiranisingh/center.html
%
% Input:
% p1 = [x1 y1]
% p2 = [x2 y2]
% p3 = [x3 y3]
%
% Output:
% circle = [x, y, radius]
%  where (x,y) is the center of the circle

function circle = find_circle (p1, p2, p3)

x1 = p1(1);
y1 = p1(2);

x2 = p2(1);
y2 = p2(2);

x3 = p3(1);
y3 = p3(2);


N1 = det([x2^2+y2^2-(x1^2+y1^2), (y2-y1); 
      x3^2+y3^2-(x1^2+y1^2), (y3-y1)]);
N2 = det([(x2-x1), x2^2+y2^2-(x1^2+y1^2);
      (x3-x1), x3^2+y3^2-(x1^2+y1^2)]);

D = det([(x2-x1), (y2-y1);
         (x3-x1), (y3-y1)]);

x = N1/(2*D);
y = N2/(2*D);

circle = [x, y, sqrt((x-x1)^2+(y-y1)^2)];
