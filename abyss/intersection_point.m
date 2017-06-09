function point = intersection_point(A, B, C, theta)
% Returns 3D coordinates of point T on side AB such that angle TCA = theta
    BAC = atan2(norm(cross(C-A, B-A)), dot(C-A, B-A));
    CTA = pi - theta - BAC;
    AT = norm(C-A) * sin(theta) / sin(CTA);
    point = A + AT.*(B-A)/norm(B-A);
end