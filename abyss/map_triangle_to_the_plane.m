function C_xy = map_triangle_to_the_plane(A, B, C)
% Returns planar coordianates for the point C such that: 1) triangles
% with planar and with spatial coordinates of A, B and C are congruent;
% 2) C.xy is on the right of A.xy -> B.xy

    ABC = angle_between_3d_vectors(A.xyz - B.xyz, C.xyz - B.xyz);
    R = [cos(ABC), -sin(ABC); sin(ABC), cos(ABC)]; % Rotation matrix
    BA_xy_unit = (A.xy-B.xy)' ./ norm(A.xy-B.xy); % Planar unit vector in direction BA
    C_xy = B.xy' + R * BA_xy_unit .* norm(C.xyz-B.xyz);
    C_xy = C_xy';
end

function angle = angle_between_3d_vectors(a, b)
    angle = atan2(norm(cross(a, b)), dot(a, b));
end