// mesh element size
lc = 1; // one micrometer
dbl = 100; // boundary layer thickness [100mu]
jnct = 10/1000; // junction thickness
memb = 100; // membrane thickness [100mu]
Point(1) = {0.0, 0.0, 0.0, lc/10};
Point(2) = {dbl, 0.0, 0.0, lc/10};
Point(3) = {dbl+memb, 0.0, 0.0, lc/10};
Point(4) = {dbl+memb+jnct, 0.0, 0.0, lc/10};
Point(5) = {dbl+memb+memb+jnct, 0.0, 0.0, lc/10};
Point(6) = {dbl+memb+memb+dbl, 0.0, 0.0, lc/10};
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};

Field[1] = BoundaryLayer;
Field[1].NodesList = {2, 3, 4, 5};
Field[1].hfar = lc/10;
Field[1].hwall_n = 0.0002;
Field[1].hwall_t = 0.0002;
Field[1].ratio = 1.01;

Background Field = 1;
