// mesh element size
lc = 1; // one micrometer
dbl = 100; // boundary layer thickness [1mu]
memb = 100; // membrane thickness [1mu]
poly = 0.10;
denom = 10000;
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {dbl, 0.0, 0.0, lc};
Point(3) = {dbl+memb, 0.0, 0.0, lc};
Point(4) = {dbl+memb+poly, 0.0, 0.0, lc};
Point(5) = {dbl+memb+dbl, 0.0, 0.0, lc};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
//Line(5) = {5,6};
//Line(6) = {6,7};
//Line(7) = {7,8};

Field[1] = BoundaryLayer;
Field[1].NodesList = {2, 3, 4, 5};
Field[1].hfar = lc;
Field[1].hwall_n = lc/denom;
Field[1].hwall_t = lc/denom;
Field[1].ratio = 1.01;

Background Field = 1;
