function bb = minBoundingBox(X)
% compute the minimum bounding box of a set of 2D points
%   Use:   boundingBox = minBoundingBox(point_matrix)
%
% Input:  2xn matrix containing the [x,y] coordinates of n points
%         *** there must be at least 3 points which are not collinear
% output: 2x4 matrix containing the coordinates of the bounding box corners
%
% Example : generate a random set of point in a randomly rotated rectangle
%     n = 50000;
%     t = pi*rand(1);
%     X = [cos(t) -sin(t) ; sin(t) cos(t)]*[7 0; 0 2]*rand(2,n);
%     X = [X  20*(rand(2,1)-0.5)];  % add an outlier
% 
%     tic
%     c = minBoundingBox(X);
%     toc
% 
%     figure(42);
%     hold off,  plot(X(1,:),X(2,:),'.')
%     hold on,   plot(c(1,[1:end 1]),c(2,[1:end 1]),'r')
%     axis equal

%{
Copyright (c) 2011, julien diener
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
%}





% compute the convex hull (CH is a 2*k matrix subset of X)
k = convhull(X(1,:),X(2,:));
CH = X(:,k);

% compute the angle to test, which are the angle of the CH edges as:
%   "one side of the bounding box contains an edge of the convex hull"
E = diff(CH,1,2);           % CH edges
T = atan2(E(2,:),E(1,:));   % angle of CH edges (used for rotation)
T = unique(mod(T,pi/2));    % reduced to the unique set of first quadrant angles

% create rotation matrix which contains
% the 2x2 rotation matrices for *all* angles in T
% R is a 2n*2 matrix
R = cos( reshape(repmat(T,2,2),2*length(T),2) ... % duplicate angles in T
       + repmat([0 -pi ; pi 0]/2,length(T),1));   % shift angle to convert sine in cosine

% rotate CH by all angles
RCH = R*CH;

% compute border size  [w1;h1;w2;h2;....;wn;hn]
% and area of bounding box for all possible edges
bsize = max(RCH,[],2) - min(RCH,[],2);
area  = prod(reshape(bsize,2,length(bsize)/2));

% find minimal area, thus the index of the angle in T 
[a,i] = min(area);

% compute the bound (min and max) on the rotated frame
Rf    = R(2*i+[-1 0],:);   % rotated frame
bound = Rf * CH;           % project CH on the rotated frame
bmin  = min(bound,[],2);
bmax  = max(bound,[],2);

% compute the corner of the bounding box
Rf = Rf';
bb(:,4) = bmax(1)*Rf(:,1) + bmin(2)*Rf(:,2);
bb(:,1) = bmin(1)*Rf(:,1) + bmin(2)*Rf(:,2);
bb(:,2) = bmin(1)*Rf(:,1) + bmax(2)*Rf(:,2);
bb(:,3) = bmax(1)*Rf(:,1) + bmax(2)*Rf(:,2);

