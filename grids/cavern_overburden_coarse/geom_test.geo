// Note: cavern height is 224.82 meters.

coarse_size_1 = 180;
coarse_size_2 = 25;
fine_size = 4.5;

Lz = 660.0;
Ly = 2000.0;
Lx = Ly;

ovb_thickness = 400;
salt_thickness = 1000;
hanging_wall = 400;

depth = ovb_thickness + hanging_wall;
depth_aux = 430 + depth;

// Cavern profile
Point(0) = {0.000000, 0.000000, 430.000000 - depth_aux, fine_size};
Point(1) = {0.000000, 0.000000, 205.189718 - depth_aux, fine_size};
Point(2) = {0.000000, 0.000000, 221.346389 - depth_aux, coarse_size_1};
Point(3) = {0.000000, 0.000000, 424.920441 - depth_aux, coarse_size_1};
Point(4) = {0.000000, 0.000000, 235.887393 - depth_aux, coarse_size_1};
Point(5) = {0.000000, 0.000000, 393.414933 - depth_aux, coarse_size_1};
Point(6) = {0.000000, 0.000000, 311.823745 - depth_aux, coarse_size_1};
Point(7) = {0.000000, 0.000000, 301.321909 - depth_aux, coarse_size_1};
Point(8) = {0.000000, 0.000000, 380.489596 - depth_aux, coarse_size_1};
Point(9) = {0.000000, 0.000000, 344.944920 - depth_aux, coarse_size_1};
Point(10) = {0.000000, 0.000000, 412.802938 - depth_aux, coarse_size_1};
Point(11) = {0.000000, 0.000000, 227.809058 - depth_aux, coarse_size_1};
Point(12) = {0.000000, 0.000000, 402.301102 - depth_aux, coarse_size_1};
Point(13) = {0.000000, 0.000000, 331.211750 - depth_aux, coarse_size_1};
Point(14) = {0.000000, 0.000000, 251.236230 - depth_aux, coarse_size_1};
Point(15) = {0.000000, 0.000000, 356.254590 - depth_aux, coarse_size_1};
Point(16) = {0.000000, 0.000000, 285.165239 - depth_aux, coarse_size_1};
Point(17) = {0.000000, 0.000000, 267.392901 - depth_aux, coarse_size_1};
Point(28) = {45.000000, 0.000000, 344.944920 - depth_aux, fine_size};
Point(29) = {68.597561, 0.000000, 251.236230 - depth_aux, fine_size};
Point(30) = {47.195122, 0.000000, 331.211750 - depth_aux, fine_size};
Point(31) = {72.987805, 0.000000, 285.165239 - depth_aux, fine_size};
Point(32) = {47.743902, 0.000000, 235.887393 - depth_aux, fine_size};
Point(33) = {47.743902, 0.000000, 380.489596 - depth_aux, fine_size};
Point(34) = {0.000000, 48.292683, 356.254590 - depth_aux, fine_size};
Point(35) = {48.292683, 0.000000, 356.254590 - depth_aux, fine_size};
Point(36) = {0.000000, 42.804878, 393.414933 - depth_aux, fine_size};
Point(37) = {0.000000, 68.597561, 251.236230 - depth_aux, fine_size};
Point(38) = {0.000000, 56.524390, 311.823745 - depth_aux, fine_size};
Point(39) = {0.000000, 74.634146, 267.392901 - depth_aux, fine_size};
Point(40) = {0.000000, 45.000000, 344.944920 - depth_aux, fine_size};
Point(41) = {0.000000, 72.987805, 285.165239 - depth_aux, fine_size};
Point(42) = {56.524390, 0.000000, 311.823745 - depth_aux, fine_size};
Point(43) = {7.134146, 0.000000, 221.346389 - depth_aux, fine_size};
Point(44) = {34.573171, 0.000000, 402.301102 - depth_aux, fine_size};
Point(45) = {0.000000, 19.756098, 412.802938 - depth_aux, fine_size};
Point(46) = {0.000000, 47.195122, 331.211750 - depth_aux, fine_size};
Point(47) = {57.621951, 0.000000, 301.321909 - depth_aux, fine_size};
Point(48) = {74.634146, 0.000000, 267.392901 - depth_aux, fine_size};
Point(49) = {0.000000, 7.134146, 221.346389 - depth_aux, fine_size};
Point(50) = {0.000000, 10.426829, 424.920441 - depth_aux, fine_size};
Point(51) = {19.756098, 0.000000, 412.802938 - depth_aux, fine_size};
Point(52) = {0.000000, 21.951220, 227.809058 - depth_aux, fine_size};
Point(53) = {0.000000, 47.743902, 235.887393 - depth_aux, fine_size};
Point(54) = {0.000000, 47.743902, 380.489596 - depth_aux, fine_size};
Point(55) = {10.426829, 0.000000, 424.920441 - depth_aux, fine_size};
Point(56) = {42.804878, 0.000000, 393.414933 - depth_aux, fine_size};
Point(57) = {21.951220, 0.000000, 227.809058 - depth_aux, fine_size};
Point(58) = {0.000000, 57.621951, 301.321909 - depth_aux, fine_size};
Point(59) = {0.000000, 34.573171, 402.301102 - depth_aux, fine_size};

// Overburden layer
Point(60) = {0.0, 0.0, 0.0, coarse_size_2};
Point(61) = {Lx, 0.0, 0.0, coarse_size_1};
Point(62) = {Lx, Ly, 0.0, coarse_size_1};
Point(63) = {0.0, Ly, 0.0, coarse_size_1};
Point(64) = {0.0, 0.0, -ovb_thickness, coarse_size_2};
Point(65) = {Lx, 0.0, -ovb_thickness, coarse_size_1};
Point(66) = {Lx, Ly, -ovb_thickness, coarse_size_1};
Point(67) = {0.0, Ly, -ovb_thickness, coarse_size_1};


// Salt layer
Point(68) = {0.0, 0.0, -ovb_thickness-salt_thickness, coarse_size_2};
Point(69) = {Lx, 0.0, -ovb_thickness-salt_thickness, coarse_size_1};
Point(70) = {Lx, Ly, -ovb_thickness-salt_thickness, coarse_size_1};
Point(71) = {0.0, Ly, -ovb_thickness-salt_thickness, coarse_size_1};

Spline(1) = {1, 43, 57, 32, 29, 48, 31, 47, 42};
Spline(2) = {1, 49, 52, 53, 37, 39, 41, 58, 38};

Spline(3) = {42, 30, 28, 35, 33, 56, 44, 51, 55, 0};
Spline(4) = {38, 46, 40, 34, 54, 36, 59, 45, 50, 0};

Line(3) = {1, 68};
Line(4) = {68, 69};
Line(5) = {69, 70};
Line(6) = {70, 71};
Line(7) = {71, 67};
Line(8) = {67, 64};
Line(9) = {64, 65};
Line(10) = {65, 66};
Line(11) = {66, 67};
Line(12) = {67, 63};
Line(13) = {63, 60};
Line(14) = {60, 61};
Line(15) = {61, 62};
Line(16) = {62, 63};
Line(17) = {62, 66};
Line(18) = {66, 70};
Line(19) = {69, 65};
Line(20) = {65, 61};
Line(21) = {71, 68};
Line(22) = {60, 64};
Line(23) = {64, 0};
Line(25) = {1, 68};
Line(26) = {68, 69};

Circle(24) = {42, 6, 38};
Curve Loop(1) = {1, 24, -2};
Surface(1) = {1};
Curve Loop(2) = {24, 4, -3};
Surface(2) = {2};
//+
Curve Loop(3) = {26, 19, -9, 23, -3, -1, 25};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {25, -21, 7, 8, 23, -4, -2};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {21, 26, 5, 6};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {9, 20, -14, 22};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {22, -8, 12, 13};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {11, 12, -16, 17};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {17, -10, 20, 15};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {15, 16, 13, 14};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {11, 8, 9, 10};
//+
Plane Surface(11) = {11};
//+
Curve Loop(12) = {10, 18, -5, 19};
//+
Plane Surface(12) = {12};
//+
Curve Loop(13) = {18, 6, 7, -11};
//+
Plane Surface(13) = {13};
//+
Surface Loop(1) = {8, 7, 6, 9, 10, 11};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {5, 4, 3, 12, 13, 2, 1, 11};
//+
Volume(2) = {2};
