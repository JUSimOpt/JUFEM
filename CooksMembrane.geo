/* Script to generate Cooks Membrane */

Point(1) = {0, 0, 0, 1.0};
Point(2) = {48, 44, 0, 1.0};
Point(3) = {48, 60, 0, 1.0};
Point(4) = {0, 44, 0, 1.0};
Point(5) = {0, 0, 16, 1.0};
Point(6) = {48, 44, 16, 1.0};
Point(7) = {48, 60, 16, 1.0};
Point(8) = {0, 44, 16, 1.0};
Line(1) = {8, 5};
Line(2) = {5, 6};
Line(3) = {6, 7};
Line(4) = {7, 8};
Line(5) = {5, 1};
Line(6) = {1, 2};
Line(7) = {2, 6};
Line(8) = {2, 3};
Line(9) = {3, 7};
Line(10) = {3, 4};
Line(11) = {4, 8};
Line(12) = {4, 1};


// Lines with dimensions type, 4 points
Transfinite Line{3,8} = 4 Using Progression 1;

// Lines with dimensions type, 2 points
Transfinite Line{5,7,9,11} = 3 Using Progression 1;

// Lines with dimensions type, 4 points
Transfinite Line{1,12} = 4 Using Progression 1;

// Lines with dimensions type, 6 points
Transfinite Line{2,4,6,10} = 6 Using Progression 1;

Line Loop(1) = {-1,-11,12,-5};
Surface(1) = {1};
Transfinite Surface{1};
Recombine Surface{1}; // Recombine the triangles into quads

Line Loop(2) = {3,-9,-8,7};
Surface(2) = {2};
Transfinite Surface{2};
Recombine Surface{2}; // Recombine the triangles into quads

Line Loop(3) = {-4,-9,10,11};
Surface(3) = {3};
Transfinite Surface{3};
Recombine Surface{3}; // Recombine the triangles into quads

Line Loop(4) = {2,-7,-6,-5};
Surface(4) = {4};
Transfinite Surface{4};
Recombine Surface{4}; // Recombine the triangles into quads

Line Loop(5) = {6,8,10,12};
Surface(5) = {5};
Transfinite Surface{5};
Recombine Surface{5}; // Recombine the triangles into quads

Line Loop(6) = {2,3,4,1};
Surface(6) = {6};
Transfinite Surface{6};
Recombine Surface{6}; // Recombine the triangles into quads

Surface Loop(1) = {3,5,4,6,1,2};
Volume(1) = {1};
Transfinite Volume{1};
Recombine Surface{6}; // Recombine the triangles into quads

//Creates 3D Mesh
Mesh 3;