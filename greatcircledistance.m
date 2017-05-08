function arclen = greatcircledistance (theta_1, phi_1, theta_2, phi_2, radius)
% From http://math.stackexchange.com/questions/231221/great-arc-distance-between-two-points-on-a-unit-sphere
% dsigma = acos( cos(theta_1) * cos(theta_2) + sin(theta_1) * sin(theta_2) * cos(phi_1 - phi_2) );
% arclen = radius*dsigma;

[x,y,z] = sph2cart(theta_1, phi_1, radius);
xyz_1 = [x,y,z];
[x,y,z] = sph2cart(theta_2, phi_2, radius);
xyz_2 = [x,y,z];
aa = cross(xyz_1, xyz_2);
bb = dot(xyz_1', xyz_2');
aa = arrayfun(@(n) norm(aa(n, :)), 1:size(aa, 1));
arclen = radius * atan2(aa, bb)';



