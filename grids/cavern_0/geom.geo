gridsize = 45;

Lz = 660;
Ly = 450;

h = 200;
H = 170;
R = 45;


Point(1) = {0.0, 0.0, 0.0, gridsize};
Point(2) = {Ly, 0.0, 0.0, gridsize};
Point(3) = {Ly, 0.0, -Lz, gridsize};
Point(4) = {0.0, 0.0, -Lz, gridsize};

Point(5) = {0.0, 0.0, -h, gridsize/10};
Point(6) = {0.0, 0.0, -h-R, gridsize};
Point(7) = {0.0, 0.0, -h-R-H, gridsize};
Point(8) = {0.0, 0.0, -h-R-H-R, gridsize/10};

Point(9) = {R, 0.0, -h-R, gridsize/10};
Point(10) = {R, 0.0, -h-R-H, gridsize/10};

Point(11) = {0, R, -h-R-H, gridsize/10};
Point(12) = {0, R, -h-R, gridsize/10};
Point(13) = {0, Ly, 0, gridsize};
Point(14) = {0, Ly, -Lz, gridsize};



Line(1) = {8, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line(5) = {1, 5};
Circle(6) = {5, 6, 9};
Circle(7) = {8, 7, 10};
Line(8) = {9, 10};

Line(9) = {4, 14};
Line(10) = {14, 13};
Line(11) = {13, 1};
Circle(12) = {14, 4, 3};
Circle(13) = {13, 1, 2};
Circle(14) = {8, 7, 11};
Circle(15) = {11, 7, 10};
Circle(16) = {5, 6, 12};
Circle(17) = {12, 6, 9};
Line(18) = {12, 11};


Curve Loop(1) = {12, 3, -13, -10};
Surface(1) = {-1};
Curve Loop(2) = {4, -11, 13};
Plane Surface(2) = {-2};
Curve Loop(3) = {5, 6, 8, -7, 1, 2, 3, 4};
Plane Surface(3) = {3};
Curve Loop(4) = {9, 10, 11, 5, 16, 18, -14, 1};
Plane Surface(4) = {-4};
Curve Loop(5) = {9, 12, -2};
Plane Surface(5) = {5};
Curve Loop(6) = {7, -15, -14};
Surface(6) = {6};
Curve Loop(7) = {15, -8, -17, 18};
Surface(7) = {7};
Curve Loop(8) = {17, -6, 16};
Surface(8) = {8};

Surface Loop(1) = {3, 4, 6, 7, 8, 2, 1, 5};
Volume(1) = {1};


Physical Surface("BOTTOM", 19) = {5};
Physical Surface("TOP", 20) = {2};
Physical Surface("OUTER", 21) = {1};
Physical Surface("SIDE_X", 22) = {4};
Physical Surface("SIDE_Y", 23) = {3};
Physical Surface("WALL", 24) = {6, 7, 8};
Physical Volume("BODY", 25) = {1};
