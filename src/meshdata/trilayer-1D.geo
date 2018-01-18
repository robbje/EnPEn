// mesh element size
lc = 1; // one micrometer
thickness = 100.0; // boundary layer thickness [1mu]
memb_thickness = 100.0; // membrane thickness [1mu]
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {thickness, 0.0, 0.0, lc};
Point(3) = {thickness+memb_thickness, 0.0, 0.0, lc};
Point(4) = {thickness+memb_thickness+thickness, 0.0, 0.0, lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};

Field[1] = BoundaryLayer;
Field[1].NodesList = {2, 3};
Field[1].hfar = lc;
Field[1].hwall_n = lc/100000;
Field[1].hwall_t = lc/100000;
Field[1].ratio = 1.01;

Background Field = 1;
