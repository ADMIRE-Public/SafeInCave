r0 = 0.001;
r1 = 0.3;

mm = 1/1000;

Lx = 48.0*mm;
Ly1 = 44.0*mm;
Ly2 = 16.0*mm;
Lz = 1*mm;


Point(1) = {0.0, 0.0, 0.0, r0};
Point(2) = {0.0, Ly1, 0.0, r0};
Point(3) = {Lx, Ly1, 0.0, r0};
Point(4) = {Lx, Ly1+Ly2, 0.0, r0};

Point(5) = {0.0, 0.0, Lz, r0};
Point(6) = {0.0, Ly1, Lz, r0};
Point(7) = {Lx, Ly1, Lz, r0};
Point(8) = {Lx, Ly1+Ly2, Lz, r0};

Line(1) = {1, 5};
Line(2) = {5, 7};
Line(3) = {7, 3};
Line(4) = {3, 1};
Line(5) = {1, 2};
Line(6) = {2, 6};
Line(7) = {6, 8};
Line(8) = {8, 4};
Line(9) = {4, 3};
Line(10) = {6, 5};
Line(11) = {2, 4};
Line(12) = {8, 7};

Curve Loop(1) = {4, 5, 11, 9};
Plane Surface(1) = {1};
Curve Loop(2) = {10, 2, -12, -7};
Plane Surface(2) = {2};
Curve Loop(3) = {7, 8, -11, 6};
Plane Surface(3) = {3};
Curve Loop(4) = {10, -1, 5, 6};
Plane Surface(4) = {-4};
Curve Loop(5) = {1, 2, 3, 4};
Plane Surface(5) = {-5};
Curve Loop(6) = {9, -3, -12, 8};
Plane Surface(6) = {-6};
Surface Loop(1) = {2, 4, 5, 6, 1, 3};
Volume(1) = {1};

Physical Surface("TOP", 13) = {2};
Physical Surface("BOTTOM", 14) = {1};
Physical Surface("BOTTOM", 14) += {1};
Physical Surface("WEST", 15) = {4};
Physical Surface("EAST", 16) = {6};
Physical Surface("SOUTH", 17) = {5};
Physical Surface("NORTH", 18) = {3};
//+
Physical Volume("BODY", 19) = {1};
