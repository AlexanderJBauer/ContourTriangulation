// Alexander Bauer
// CS 130B Programming Assignment 1
// Last Edited: 4/24/2017

#include <string>   /* std::string */
#include <math.h> /* sqrt */
#include <vector> /* std::vector */
#include <algorithm> /* std::min */
#include <iostream> /* std::cout, std::cin */
#include <limits> /* std::numeric_limits */
#include <iomanip> /* std::setw */

int numComparisons = 0; //Global variable for counting comparisons

// Point struct to make things easier
struct Point{

    double x;
    double y;
    double z;


    //CONSTRUCTOR
    Point( double xVal = 0, double yVal = 0, double zVal = 0 )
    : x{ xVal }, y{ yVal }, z{ zVal } { }

};

// DISTANCE FUNCTION
double distanceBetween( Point p1, Point p2 ){
    return sqrt( ( p1.x - p2.x ) * ( p1.x - p2.x )
                +( p1.y - p2.y ) * ( p1.y - p2.y )
                +( p1.z - p2.z ) * ( p1.z - p2.z ) );
}

// AREA FUNCTION
double areaOf( Point p1, Point p2, Point p3 ){
    double a = distanceBetween( p1, p2 );
    double b = distanceBetween( p2, p3 );
    double c = distanceBetween( p3, p1 );
    double s = ( a + b + c )/2;

    return sqrt( s * ( s-a ) * ( s-b ) * ( s-c ) );
}

double MinimumCostPath( std::vector<Point> pointsOfP, std::vector<Point> pointsOfQ, int s ){

    int q = pointsOfQ.size();
    int p = pointsOfP.size();
    int numRows =  q + 1;
    int numColumns = p + 1;

    std::vector< std::vector<double> > costMatrix( numRows, std::vector<double>(numColumns) );

    // Calculate top row
    for( int i = 1; i < numColumns; i++ ){
        costMatrix[0][i] = costMatrix[0][i-1] + areaOf( pointsOfQ[0], pointsOfP[(i-1+s)%p], pointsOfP[(i+s)%p] );

    }

    // Calculate leftmost column
    for( int i = 1; i < numRows; i++ ){
        costMatrix[i][0] = costMatrix[i-1][0] + areaOf( pointsOfP[0+s], pointsOfQ[i-1], pointsOfQ[i%q] );
    }

    // Rest of Matrix
    for( int row = 1; row < numRows; row++ ){
        for( int col = 1; col < numColumns; col++ ){
            costMatrix[row][col] = std::min(
              (costMatrix[row][col-1] + areaOf( pointsOfQ[row%q], pointsOfP[(col-1+s)%p], pointsOfP[(col+s)%p] )),
              (costMatrix[row-1][col] + areaOf( pointsOfP[(col+s)%p], pointsOfQ[row-1], pointsOfQ[row%q] ))
            );
        }
    }

    return costMatrix[q][p];
}

std::vector< std::vector<int> > CreatePath( std::vector<Point> pointsOfP, std::vector<Point> pointsOfQ, int s ){
        int q = pointsOfQ.size();
        int p = pointsOfP.size();
        int numRows =  q + 1;
        int numColumns = p + 1;

        std::vector< std::vector<double> > costMatrix( numRows, std::vector<double>(numColumns) );

        // Calculate top row
        for( int i = 1; i < numColumns; i++ ){
            costMatrix[0][i] = costMatrix[0][i-1] + areaOf( pointsOfQ[0], pointsOfP[(i-1+s)%p], pointsOfP[(i+s)%p] );

        }

        // Calculate leftmost column
        for( int i = 1; i < numRows; i++ ){
            costMatrix[i][0] = costMatrix[i-1][0] + areaOf( pointsOfP[0+s], pointsOfQ[i-1], pointsOfQ[i%q] );
        }

        // Rest of Matrix
        for( int row = 1; row < numRows; row++ ){
            for( int col = 1; col < numColumns; col++ ){
                costMatrix[row][col] = std::min(
                  (costMatrix[row][col-1] + areaOf( pointsOfQ[row%q], pointsOfP[(col-1+s)%p], pointsOfP[(col+s)%p] )),
                  (costMatrix[row-1][col] + areaOf( pointsOfP[(col+s)%p], pointsOfQ[row-1], pointsOfQ[row%q] ))
                );
            }
        }

        std::vector< std::vector<int> > triangles( p+q, std::vector<int>(3) );

        int row = q;
        int col = p;

        for( int i = 0; i < p + q; i++ ){
            if( row == 0 ){
                triangles[i][0] = row % q + p;
                triangles[i][1] = (col-1+s)%p;
                triangles[i][2] = (col+s)%p;
                col = col-1;
            }
            else if( col == 0 ){
                triangles[i][0] = (col+s) % p;
                triangles[i][1] = (row-1) + p;
                triangles[i][2] = row % q + p;
                row = row-1;
            }
            else if( costMatrix[row-1][col] < costMatrix[row][col-1] ){
                triangles[i][0] = (col+s) % p;
                triangles[i][1] = (row-1) + p;
                triangles[i][2] = row % q + p;
                row = row-1;
            }
            else{
                triangles[i][0] = row % q + p;
                triangles[i][1] = (col-1+s)%p;
                triangles[i][2] = (col+s)%p;
                col = col-1;
            }
        }
        return triangles;
}

int main() {

    int numPointsP = 0;
    int numPointsQ = 0;
    double x;
    double y;
    double z;

    std::cin >> numPointsP;
    std::cin >> numPointsQ;

    std::vector<Point> allPoints;

    std::vector<Point> pointsOfP;
    for( int i = 0; i < numPointsP; i++){
        std::cin >> x;
        std::cin >> y;
        std::cin >> z;
        pointsOfP.push_back( Point( x, y, z ) );
        allPoints.push_back( Point( x, y, z ) );
    }

    std::vector<Point> pointsOfQ;
    for( int i = 0; i < numPointsQ; i++){
        std::cin >> x;
        std::cin >> y;
        std::cin >> z;
        pointsOfQ.push_back( Point( x, y, z ) );
        allPoints.push_back( Point( x, y, z ) );
    }

        std::cout.precision(4);
        double knownCost = 0;
        int a = 0;
        int b = 0;
        int c = 0;
        for( int i = 0; i < numPointsP + numPointsQ; i++ ){
            std::cin >> a;
            std::cin >> b;
            std::cin >> c;

            knownCost = knownCost + areaOf( allPoints[a-1],allPoints[b-1],allPoints[c-1] );
        }

    double minCost = std::numeric_limits<double>::max();
    double newCost = 0;
    int bestIndex = 0;
    for( int i = 0; i < numPointsP; i++ ){
        newCost = MinimumCostPath( pointsOfP, pointsOfQ, i );
        if( newCost < minCost ){
            minCost = newCost;
            bestIndex = i;
        }
    }

    std::cout << minCost << std::endl << bestIndex << std:: endl;

    std::cout << knownCost << std::endl;

    std::vector< std::vector<int> > triangles = CreatePath( pointsOfP, pointsOfQ, bestIndex );

    for( int i = 0; i < allPoints.size(); i ++ ){
        std::cout << std::setw(3) << triangles[i][0]+1 << std::setw(6) << triangles[i][1]+1 << std::setw(6) << triangles[i][2]+1 << std::endl;
    }

    return 0;
}
