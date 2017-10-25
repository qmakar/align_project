#include "align.h"
#include <string>
#include "matrix.h"
#include "io.h"

using std::string;
using std::cout;
using std::endl;
using std::tuple;
using std::get;
using std::tie;
using std::make_tuple;


class Filter
{
    // Radius of neighbourhoud, which is passed to that operator
public:
    Matrix<double> kernel;
    int radius;
    uint rows, cols;


    Filter(Matrix<double>& matrix, int rad, int r, int c): kernel(matrix), radius(rad), rows(r), cols(c) {}
    tuple<uint, uint, uint> operator () (const Image &m) const {
        // uint size = 2 * radius + 1;
        uint r, g, b, sum_r = 0, sum_g = 0, sum_b = 0;
        for (uint i = 0; i < rows; ++i) {
            for (uint j = 0; j < cols; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = m(i, j);
                sum_r += r * kernel(i,j);
                sum_g += g * kernel(i,j);
                sum_b += b * kernel(i,j);
            }
        }
        return make_tuple(sum_r, sum_g, sum_b);
    }
};




Image align(Image srcImage, bool isPostprocessing, std::string postprocessingType, double fraction, bool isMirror, 
            bool isInterp, bool isSubpixel, double subScale)
{
    return srcImage;
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
    return src_image;
}

Image gray_world(Image src_image) {
    return src_image;
}

Image resize(Image src_image, double scale) {
    return src_image;
}

Image custom(Image src_image, Matrix<double> kernel) {
    // Function custom is useful for making concrete linear filtrations
    // like gaussian or sobel. So, we assume that you implement custom
    // and then implement other filtrations using this function.
    // sobel_x and sobel_y are given as an example.
    
    Filter f(kernel, (kernel.n_rows - 1) / 2, kernel.n_rows, kernel.n_cols);
    return src_image.unary_map(f);
}


Image autocontrast(Image src_image, double fraction) {
    return src_image;
}


void gauss_matrix(Matrix<double>& GM, double sigma, int radius){
    int size = 2 * radius + 1;
    
    // Compute matrix by normal distribution
    int x, y;
    double si = 2.0 * sigma * sigma;
    double norm = 0;
    for (int i = 0; i < size; ++i){
        for (int j = 0; j < size; ++j){
            y = i - radius;
            x = j - radius;
            GM(i, j) = exp( - (x * x + y * y) / si );
            norm += GM(i ,j);
        }
    }

    // Normalize matrix to get sum of elements == 1
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            GM(i, j) /= norm;
        }
    }
}


Image gaussian(Image src_image, double sigma, int radius) {

    Matrix<double> GM(2 * radius + 1, 2 * radius + 1);
    gauss_matrix(GM, sigma, radius);

    uint size = 2 * radius + 1;
    double r_convolution = 0, g_convolution = 0, b_convolution = 0;
    uint r, g, b;
    
    Image tmp = src_image;
    // Matrix<double> sub = src_image.submatrix(radius, radius, src_image.n_rows - radius, src_image.n_rows - radius);
    // Image dst_image = Mirror( src_image, GM.n_rows, GM.n_cols );
    for( uint i = radius; i < src_image.n_rows - radius - 1; i++ ){
        for ( uint j = radius; j < src_image.n_cols - radius - 1; j++){
            r_convolution = g_convolution = b_convolution = 0;
            for( uint x = 0; x < size; x++){
                for( uint y = 0; y < size; y++){
                    tie(r, g, b) = tmp(i - radius + x, j - radius + y );
                    r_convolution += GM( x, y ) *  r;
                    g_convolution += GM( x, y ) *  g;
                    b_convolution += GM( x, y ) *  b;
                }
            }
            // to0or255( r_convolution, g_convolution, b_convolution );
            tmp(i, j) = make_tuple( r_convolution, g_convolution, b_convolution );
        }
    }
    return tmp.submatrix(radius, radius, src_image.n_rows - 2 * radius - 1, src_image.n_cols - 2 * radius - 1);;




    // // Porque no funciona???????????????????????????????
    // // Create Gauss matrix
    
    // // // Set filter to convolution
    // Filter f(GM, radius, GM.n_rows, GM.n_cols);
    
    // // Apply this filter to image
    // Image img = src_image.unary_map(f);
    // return img;
}


void gauss_column(Matrix<double>& GM, double sigma, int radius){
    
    // Compute matrix by normal distribution
    int y;
    double si = 2.0 * sigma * sigma;
    double norm = 0;
    for (uint i = 0; i < GM.n_rows; ++i){
        for(uint j = 0; j < GM.n_cols; ++j) {
            y = i + j - radius;                     // because i or j == 0
            GM(i, j) = exp( - (y * y) / si );
            norm += GM(i ,j);
        }
    }

    // Normalize matrix to get sum of elements == 1
    for (uint i = 0; i < GM.n_rows; ++i){
        for(uint j = 0; j < GM.n_cols; ++j) {
            GM(i, j) /= norm;
        }
    }   
}

Image gaussian_separable(Image src_image, double sigma, int radius) {

    // Create 2 Gauss lists
    Matrix<double> GC(2 * radius + 1, 1), GC2(1, 2 * radius + 1);
    gauss_column(GC, sigma, radius);
    gauss_column(GC2, sigma, radius);

    // Set filter to convolution
    Filter f(GC, radius, GC.n_rows, 1);
    Filter g(GC2, radius, 1, GC2.n_cols);
    
    // Apply this filter to image
    Image img = src_image.unary_map(f);
    img = img.unary_map(g);

    return img;
}


Matrix<tuple<double, double>> gradient(Image img1, Image img2){
    Matrix<tuple<double, double>> grad(img1.n_rows, img1.n_cols);
    uint pixel1, pixel2;
    for(unsigned i = 0; i < img1.n_rows; ++i) {
        for(unsigned j = 0; j < img1.n_cols; ++j) {
            pixel1 = get<0>(img1(i,j));
            pixel2 = get<0>(img2(i,j));
            grad(i,j) = make_tuple(atan2(pixel1, pixel2), sqrt(pixel1 * pixel1 + pixel2 * pixel2));
        }
    }
    return grad;
}

Image suppression(Matrix<tuple<double, double>> grad){
    Image tmp( grad.n_rows, grad.n_cols );
    double teta, mod;
    double pi8 = 3.14 / 8;
    int i1 = 0, i2 = 0, j1 = 0, j2 = 0;
    for ( uint i = 1; i < grad.n_rows - 1; i++ ){
        for ( uint j = 1; j < grad.n_cols - 1; j++ ){
            teta = get<0>( grad(i,j));
            if ( ( teta >= -pi8 || teta <= pi8 ) || ( teta >= 7*pi8 || teta <= -7*pi8 ) ) { i1 = i; j1 = j+1; i2 = i; j2 = j-1;};
            if ( ( teta >= pi8 && teta <= 3*pi8 ) || ( teta >= -7*pi8 && teta <= -5*pi8 ) ) { i1 = i-1; j1 = j+1; i2 = i +1; j2 = j-1;};
            if ( ( teta >= 3*pi8 && teta <= 5*pi8 ) || ( teta >= -5*pi8 && teta <= -3*pi8 ) ) { i1 = i-1; j1 = j; i2 = i+1; j2 = j;};
            if ( ( teta >= 5*pi8 && teta <= 7*pi8 ) || ( teta >= -3*pi8 && teta <= -pi8 ) ) { i1 = i-1; j1 = j-1; i2 = i+1; j2 = j+1;};
            mod = get<1>(grad(i, j));
            if ( ( mod <= get<1>(grad(i1, j1)) ) || ( mod <= get<1>(grad(i2, j2)) ) ){
                mod = 0;
            }
            if ( mod > 255 ) mod = 255;
            tmp( i, j ) = make_tuple( mod, mod, mod );
        }
    }
    return tmp;
}

Image cut_to_bounds( Image src_image, int threshold1, int threshold2){
    int r;
    for ( uint i = 0; i < src_image.n_rows; i++ ){
        for( uint j = 0; j < src_image.n_cols; j++ ){
            tie( r, r, r ) = src_image( i, j );
            if ( r < threshold1 ){
                src_image( i, j ) = make_tuple( 0, 0, 0 );
            }
            else if ( r > threshold2 ){
                src_image( i, j ) = make_tuple( 255, 255, 255);
            }
            else {
                src_image( i, j ) = make_tuple( 127, 127, 127);
            }
        }
    }
    return src_image;
}


Image median(Image src_image, int radius) {
    return src_image;
}

Image median_linear(Image src_image, int radius) {
    return src_image;
}

Image median_const(Image src_image, int radius) {
    return src_image;
}

tuple < int, int > & Find_Set( Matrix< tuple < int, int, int>> & m, tuple < int, int > & ij){
    int color, i, j;
    int a, b;
    tie( a, b ) = ij;
    tie( color, i, j ) = m( a, b );
    if ( i != a || j != b) {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!....****
        ij = make_tuple( i, j );
        ij = Find_Set( m, ij );
    }
    return ij;
}

void union_images(  Matrix< tuple < int, int, int>> & m, int a, int b, int c, int d ){
    tuple < int, int > ij1, ij2;
    ij1 = make_tuple( a, b );
    ij2 = make_tuple( c, d );
    ij1 = Find_Set( m, ij1 );
    ij2 = Find_Set( m, ij2 );
    int i1, j1, i2, j2;
    tie( i1, j1 ) = ij1;
    tie( i2, j2 ) = ij2;
    int colX, iX, jX;
    int colY, iY, jY;
    tie ( colX, iX, jX ) = m( i1, j1 );
    tie ( colY, iY, jY ) = m( i2, j2 );
    if ( colX >= colY ){// оставить как есть!!!!!!!!!!!
        iY = iX;
        jY = jX;
    }
    else{
        iX = iY;
        jX = jY;
    }
    m( a, b ) = make_tuple( colX, iX, jX );
    m( c, d ) = make_tuple( colY, iY, jY );
}

void neighbours( Matrix< tuple < int, int, int>> & m, int i, int j, Image &src_image ){
    int color, right_color, up_color, diag_color;
    tie ( color, color, color ) = src_image( i, j );
    m( i, j ) = make_tuple( color, i, j );
    right_color = get<1> ( src_image( i, j - 1 ) );
    up_color = get<1> ( src_image( i - 1, j ) );
    diag_color = get<1> ( src_image( i - 1, j - 1 ) );
    if ( up_color != 0 && right_color != 0){
        union_images( m, i - 1, j, i, j - 1);
        if ( color > 0 ) union_images( m, i, j, i, j - 1 );
    }
    else if ( color > 0){
        if( up_color > 0 ) union_images( m, i, j, i - 1, j );
        if( right_color > 0 ) union_images( m, i, j, i, j - 1 );
        if( diag_color > 0 ) union_images( m, i, j, i - 1, j - 1 );
    }
}

Image hysteresis( Image src_image ){
    Matrix< tuple < int, int, int>> component ( src_image.n_rows, src_image.n_cols );
    int color, right_color, up_color;// diag_color;
    uint i = 0, j = 0;
    tie( color, color, color ) = src_image( 0, 0 );
    component( 0, 0 ) = make_tuple( color, 0, 0);
    for(  j = 1; j < src_image.n_cols; j++){
        tie( color, color, color ) = src_image( i, j );
        component( i, j ) = make_tuple( color, i, j);
        if ( color != 0 ){
            right_color = get< 1 >( src_image( i, j - 1 ));
            if ( right_color != 0 ){
                union_images( component, i, j - 1, i, j);
            
            }
        }
    }
    j = 0;
    for(  i = 1; i < src_image.n_rows; i++){
        tie( color, color, color ) = src_image( i, j );
        component( i, j ) = make_tuple( color, i, j);
        if ( color != 0 ){
            up_color = get< 1 >( src_image( i - 1, j ));
            if ( up_color != 0 ){
                union_images( component, i - 1, j, i, j);
            
            }
        }
    }
    i = j = 1;
    uint imax = src_image.n_rows - 1, jmax = src_image.n_cols - 1;
    int min, k = 1;
    if ( imax > jmax ) min = jmax;
    else min = imax;
    while ( k <= min ){
        // для квадратной матрицы
        j = 1;
        while ( i > j ){
            neighbours( component, i, j, src_image );
            j++;
        }
        i = 1;
        while ( i <= j ){
            neighbours( component, i, j, src_image );
            i++;
        }
        k++;
    }
    int flag = 0;
    if ( i > imax ) { j++; i = imax; flag = 1;}
    if ( flag == 1 ){
        while ( j <= jmax ){ 
            for ( i = 1; i <= imax; ++i){
                neighbours( component, i, j, src_image );
            }
            j++;
        }
    }
    else{
        while ( i <= imax && j<= jmax ){ 
            for ( j = 1; j <= jmax; ++j){
                neighbours( component, i, j, src_image );
            }
            i++;
        }
    }
    int l;
    tuple< int, int > ij;
    for ( i = 0; i < src_image.n_rows; ++i){
        for ( j = 0; j < src_image.n_cols; ++j){
            ij = make_tuple( i, j );
            ij = Find_Set( component, ij );
            k = get< 0 >( ij );
            l = get< 1 >( ij );
            tie( color, color, color ) = src_image( k, l );
            if ( color == 255 )
            src_image( i, j ) = make_tuple( 255, 255, 255 );
            else
            src_image( i, j ) = make_tuple( 0, 0, 0);
        }
    }
    return src_image;
}
Image canny(Image src_image, int threshold1, int threshold2) {
    Image img = gaussian(src_image, 1.4, 2);
    Image Ix = sobel_x(img);
    Image Iy = sobel_y(img);
    Matrix<tuple<double, double>> grad = gradient(Ix, Iy);
    img = Suppression(grad);
    img = cut_to_bounds( img, threshold1, threshold2);
    img = hysteresis(img);
    return img;
}
