coarse_size = 0.1;

Lz = 1.0;
Ly = 1.0;
Lx = 1.0;

Point(0) = {0.0, 0.0, 0.0, coarse_size};
Point(1) = {Lx, 0.0, 0.0, coarse_size};
Point(2) = {Lx, Ly, 0.0, coarse_size};
Point(3) = {0.0, Ly, 0.0, coarse_size};
Point(4) = {0.0, 0.0, Lz, coarse_size};
Point(5) = {Lx, 0.0, Lz, coarse_size};
Point(6) = {Lx, Ly, Lz, coarse_size};
Point(7) = {0.0, Ly, Lz, coarse_size};

Line(1) = {4, 5};
Line(2) = {5, 1};
Line(3) = {1, 0};
Line(4) = {0, 4};
Line(5) = {4, 7};
Line(6) = {7, 6};
Line(7) = {6, 2};
Line(8) = {2, 3};
Line(9) = {3, 7};
Line(10) = {6, 5};
Line(11) = {2, 1};
Line(12) = {3, 0};

Curve Loop(1) = {1, -10, -6, -5};
Plane Surface(1) = {1};
Curve Loop(2) = {4, 5, -9, 12};
Plane Surface(2) = {2};
Curve Loop(3) = {1, 2, 3, 4};
Plane Surface(3) = {-3};
Curve Loop(4) = {2, -11, -7, 10};
Plane Surface(4) = {4};
Curve Loop(5) = {8, 12, -3, -11};
Plane Surface(5) = {-5};
Curve Loop(6) = {6, 7, 8, 9};
Plane Surface(6) = {6};

Surface Loop(1) = {6, 1, 3, 2, 5, 4};
Volume(1) = {1};

Physical Surface("TOP", 13) = {1};
Physical Surface("BOTTOM", 14) = {5};
Physical Surface("EAST", 15) = {4};
Physical Surface("WEST", 16) = {2};
Physical Surface("SOUTH", 17) = {3};
Physical Surface("NORTH", 18) = {6};
Physical Volume("BODY", 19) = {1};
