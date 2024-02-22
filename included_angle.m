function angle = included_angle(v1, v2)
    % calculate the included angle between two vectors
    dot_prod = dot(v1, v2);
    cos_theta = dot_prod / (norm(v1) * norm(v2));
    angle = acos(cos_theta);
end
