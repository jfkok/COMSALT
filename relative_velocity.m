function v_rel = relative_velocity (v_x1, v_z1, v_x2, v_z2)

%this calculates the magnitude of the relative velocity between v_x and v_z
v_rel = sqrt((v_x1-v_x2)^2+(v_z1-v_z2)^2);