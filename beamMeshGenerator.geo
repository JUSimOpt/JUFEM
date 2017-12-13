/* Script to generate beam

Mirza Cenanovic, Nov 2017
*/

L=10; //mm
W=2;  //mm
H=2;  //mm
d = 0; //mm
//d = 0.2;
//d = 0.4;

Point(1) = {0, 0, 0};
Point(2) = {L, 0, 0};
Point(3) = {L, W, 0};
Point(4) = {0, W, 0};
Point(5) = {0, 0, H};
Point(6) = {L, 0, H};
Point(7) = {L, W, H};
Point(8) = {0, W, H};
Point(9) = {5-d, 0, 0};
Point(10) = {5, 0, H};
Point(11) = {5-d, W, H};
Point(12) = {5, W, 0};

Line(1) = {1, 9};
Line(2) = {9, 2};
Line(3) = {8, 11};
Line(4) = {11, 7};
Line(5) = {4, 12};
Line(6) = {12, 3};
Line(7) = {9, 10};
Line(8) = {10, 11};
Line(9) = {11, 12};
Line(10) = {9, 12};
Line(11) = {5, 10};
Line(12) = {10, 6};
Line(13) = {1, 5};
Line(14) = {5, 8};
Line(15) = {8, 4};
Line(16) = {4, 1};
Line(17) = {2, 6};
Line(18) = {6, 7};
Line(19) = {7, 3};
Line(20) = {3, 2};

Transfinite Line{14,16,18,20,8,10} = 3 Using Progression 1; // Create 3 points on these lines
Transfinite Line{13,15,17,19,9,7} = 3 Using Progression 1; // Create 3 points on these lines
Transfinite Line{1,2,5,6,3,4,11,12} = 6 Using Progression 1; // Create 3 points on these lines

Line Loop(1) = {14, 15, 16, 13};
Surface(1) = {1};
Transfinite Surface{1} = {1,5,8,4}; // Creates a structured mesh on this surface
Recombine Surface{1}; // Recombine the triangles into quads

Line Loop(2) = {18, 17, 20, 19};
Surface(2) = {2};
Transfinite Surface{2} = {7,6,2,3}; // Creates a structured mesh on this surface
Recombine Surface{2}; // Recombine the triangles into quads

Line Loop(3) = {8, 9, -10, 7};
Surface(3) = {3};
Transfinite Surface{3} = {10,9,12,11}; // Creates a structured mesh on this surface
Recombine Surface{3}; // Recombine the triangles into quads

Line Loop(4) = {9,-5,-15,3};
Surface(4) = {4};
Transfinite Surface{4}; // Creates a structured mesh on this surface
Recombine Surface{4}; // Recombine the triangles into quads

Line Loop(5) = {1, 10, -5, 16};
Surface(5) = {5};
Transfinite Surface{5}; // Creates a structured mesh on this surface
Recombine Surface{5}; // Recombine the triangles into quads

Line Loop(6) = {11, -7, -1, 13};
Surface(6) = {6};
Transfinite Surface{6}; // Creates a structured mesh on this surface
Recombine Surface{6}; // Recombine the triangles into quads

Line Loop(7) = {3, -8, -11, 14};
Surface(7) = {7};
Transfinite Surface{7}; // Creates a structured mesh on this surface
Recombine Surface{7}; // Recombine the triangles into quads

Line Loop(8) = {12, -17, -2, 7};
Surface(8) = {8};
Transfinite Surface{8}; // Creates a structured mesh on this surface
Recombine Surface{8}; // Recombine the triangles into quads

Line Loop(9) = {2, -20, -6, -10};
Surface(9) = {9};
Transfinite Surface{9}; // Creates a structured mesh on this surface
Recombine Surface{9}; // Recombine the triangles into quads

Line Loop(10) = {4, -18, -12, 8};
Surface(10) = {10};
Transfinite Surface{10}; // Creates a structured mesh on this surface
Recombine Surface{10}; // Recombine the triangles into quads

Line Loop(11) = {4, 19, -6, -9};
Surface(11) = {11};
Transfinite Surface{11}; // Creates a structured mesh on this surface
Recombine Surface{11}; // Recombine the triangles into quads

Surface Loop(1) = {4, 5, 6, 7, 1, 3};
Volume(1) = {1};
Transfinite Volume{1};

Surface Loop(2) = {11, 10, 2, 9, 8, 3};
Volume(2) = {2};
Transfinite Volume{2};


Mesh 3; //Creates 3D Mesh

//RefineMesh;
RefineMesh;
