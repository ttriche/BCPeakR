#ifndef _myFun_h
#define _myFun_h

#include	<iostream>
#include	<fstream>
#include        <cassert>
#include	"../TNT/tnt.h"
#include	"../JAMA_C/my_jama_cholesky.h"
#include	"../JAMA_C/my_jama_lu.h"

using namespace std;
using namespace TNT;
using namespace JAMA;

/*  E = 2.7183, Pi = 3.1416 */

/***************************************************************************
 ***        I/O file 
 ***************************************************************************/

/* This function reads the vector into a file.
 */

template <class T>
static void   ReadVector(string filename, Vector<T> &vec)
{
	ifstream	infile;
	int		i=0;
	infile.open( filename.c_str() );
	if (!infile)	{ 
	    cout << "Couldn't open file!" << endl; 
	    return;
	}
	while (infile.good())	infile >> vec[i++];
	infile.close();
}

Matrix<double>   ReadMatrix(string filename)
{
  ifstream        infile;
  infile.open( filename.c_str() );
  if (!infile) {
    cout << "Couldn't open file!" << endl;
    return  Matrix<double>();
  }

  int num_rows, num_cols;
  infile >> num_rows >> num_cols ;
  Matrix<double>  mat(num_rows, num_cols);
  if (infile.good())
    for (int i=0; i<num_rows; i++) for (int j=0; j<num_cols; j++)
      infile >> mat[i][j];

  infile.close();
  return mat;
}  

Array3D<double>  ReadArray3D(string filename)
{
  ifstream  infile;
  infile.open( filename.c_str() );
  if (!infile) {
    cout << "Couldn't open file!" << endl;
    return  Array3D<double>();
  }

  int  dim1, dim2, dim3;
  infile >> dim1 >> dim2 >> dim3;
  Array3D<double> arr3d(dim1, dim2, dim3);
  if (infile.good())
    for (int i=0; i<dim1; i++) for(int j=0; j<dim2; j++) for(int k=0; k<dim3; k++)
      infile >> arr3d[i][j][k];

  infile.close();
  return arr3d;
}





/* This function writes the vector into a file.
 */

template <class T>
static void   WriteVector(string filename, const Vector<T> &vec)
{
	ofstream	outfile;

	outfile.open( filename.c_str() );
	if (!outfile)   {
	    cout << "Couldn't open file!" << endl;
	    return;
	}
	for (int i=0; i<vec.size(); i++)    outfile << vec[i] << "\n";
	outfile.close();
}

/* This function writes a matrix into a file.
 */

template <class T>  
static void   WriteMatrix(string filename, const Matrix<T> &mat)
{
	ofstream	outfile;

	outfile.open( filename.c_str() );
	if (!outfile)	{
	    cout << "Couldn't open file!" << endl; 
	    return;
	}

	outfile << mat.num_rows() << "\t" << mat.num_cols() << "\n" << endl;

	for (int i=0; i<mat.num_rows(); i++) {
	    for (int j=0; j < mat.num_cols(); j++)
		outfile << mat[i][j] << " ";
	    outfile << "\n";
	}
	outfile.close();
}

/* This function is written specially for CPAR, it output the upper
 * triangular part of (dim1, dim2)
*/

template <class T>
static void   WriteArray3D(string filename, const Array3D<T> &arr3d)
{
	ofstream	outfile;

	outfile.open( filename.c_str() );
	if (!outfile)   {
	    cout << "Couldn't open file!" << endl;
	    return;
	}

	outfile << arr3d.dim1() << "\t" << arr3d.dim2() << "\t" << 
		arr3d.dim3() << "\n" << endl;
	for (int i=0; i<arr3d.dim1(); i++) {
	  for (int j=0; j<arr3d.dim2(); j++) {
	    for (int k=0; k<arr3d.dim3(); k++) outfile<<arr3d[i][j][k]<<" ";
	    outfile << "\n";
	  }
	  outfile << "\n" << endl;
	}
}

template <class T>
static void     Initialize(Vector<T>  &vec)
{
  T   zero=vec[0]-vec[0];
  for (int i=0; i<vec.dim(); i++)   vec[i] = zero;
}

template <class T>
static void     Initialize(Matrix<T>  &mat)
{
  T   zero=mat[0][0]-mat[0][0];
  for (int i=0; i<mat.num_rows(); i++) for(int j=0;j<mat.num_cols();j++)
    mat[i][j] = zero;
}

template <class T>
static void     Initialize(Array3D<T>  &arr)
{
  T  zero=arr[0][0][0]-arr[0][0][0];
  for (int i=0; i<arr.dim1(); i++)  for(int j=0; j<arr.dim2(); j++)
    for (int k=0; k<arr.dim3(); k++)     arr[i][j][k] = zero;
}







/**************************************************************************
 ***      Vectr/Matrix Calculations
 **************************************************************************/
/* In TNT,
    Vector: +, -, *(mul_element), dot_prod, 
    Matrix: +, -, mul_element, transpose, *(mat*mat, mat*vec)

   Here,
    Vector: mul_scalar, tensor_prod, 
    Matrix: mul_scalar, mul_vecmatvec, 


PutRowVectorInMatrix, GetRowVectorInMatrix, PutMatrixInArray3D,
GetMatrixInArray3D
*/



/* This function returns 'mat * scalar'. */
Matrix<double>  mul_scalar(const Matrix<double> &mat, double scalar)
{
  Matrix<double> res(mat.num_rows(), mat.num_cols());
  for (int i=0; i<mat.num_rows(); i++)
    for (int j=0; j<mat.num_cols(); j++)
      res[i][j] = mat[i][j] * scalar;
  return res;
}

/* This function returns 'vec * scalar'. */
Vector<double>  mul_scalar(const Vector<double> &vec, double scalar)
{
  Vector<double> res(vec.size());
  for (int i=0; i<vec.size(); i++)  res[i] = vec[i] * scalar;
  return res;
}


/* This function returns the tensor product of vecA and vecB. */
Matrix<double>  tensor_prod( const Vector<double> &vecA,
			     const Vector<double> &vecB)
{
  Matrix<double>  res(vecA.size(), vecB.size());
  for (int i=0; i<vecA.size(); i++)
    for (int j=0; j < vecB.size(); j++)
      res[i][j] = vecA[i] * vecB[j];
  return res;
}


double		mul_vecmatvec(const Vector<double> &vec, 
			const Matrix<double> &mat)
{
  assert( (vec.size()==mat.num_cols())&(mat.num_rows()==vec.size()) );
  double  res=0.0;
  for (int i=0; i<vec.size(); i++)    for (int j=0; j<vec.size(); j++)
    res += vec[i]*mat[i][j]*vec[j];
  return res;
}


/*  Matrix multiplication:  mat1 * mat2 * mat3 */
Matrix<double>   MatMulti(Matrix<double>  mat1, Matrix<double>  mat2,
			  Matrix<double>  mat3)
{
  assert( (mat1.num_cols()==mat2.num_rows()) &&
	  (mat2.num_cols()==mat3.num_rows()) );
  Matrix<double>   res(mat1.num_rows(),  mat3.num_cols());
  for (int i=0; i<mat1.num_rows(); i++)
    for (int j=0; j<mat3.num_cols(); j++) {
      res[i][j] = 0.0;
      for (int k=0; k<mat2.num_rows(); k++)
	for (int s=0; s<mat2.num_cols(); s++)
	  res[i][j] += mat1[i][k] * mat2[k][s] * mat3[s][j];
    }
  return  res;
}


/* This function returns the inverse of the matrix. If it is not square matrix, abort.
 * (Write on 05_Jul_04)
 */
Matrix<double> inv(const Matrix<double> mat)
{
  assert( mat.num_rows() == mat.num_cols() );
  Cholesky   chol(mat);
  Matrix<double>  id(mat.num_rows(), mat.num_rows());
  for (int i=0; i<mat.num_rows(); i++) for(int j=0; j<mat.num_rows(); j++)
    id[i][j] = 0.0;
  for (int i=0; i<mat.num_rows(); i++)  id[i][i] = 1.0;
  return  chol.solve(id);
}

/* This function returns the diterminant of the matrix */
double	det(const Matrix<double> mat)
{
  LU	  lu(mat);
  return  lu.det();
}

/*
double	trace(const Matrix<double> &mat)
{
  double value = 0.0;
  for (int i=0; i<mat.num_rows(); i++)  trace += mat[i][i];
  return value;
}
*/



/* These four functions put(get) rowVector (matrix) in matrix (array3d). 
   Notice the position variables are in the C++ convention. */
template <class T>
static void   PutRowVectorInMatrix(Matrix<T> &mat, int rowNum, 
				   Vector<T> &vec)
{
  assert(mat.num_cols() == vec.size());
  for (int j=0; j<mat.num_cols(); j++)
    mat[rowNum][j] = vec[j];
}

template <class T>  
Vector<T>   GetRowVectorInMatrix(Matrix<T> &mat, int rowNum)
{
  Vector<T>  vec(mat.num_cols());
  for (int i=0; i<mat.num_cols(); i++)
    vec[i] = mat[rowNum][i];
  return vec;
}

/* This function assign matrix in Array3D at index [position][*][*] */
template <class T>  
static void   PutMatrixInArray3D(Array3D<T> &arr3d, int position,
				 Matrix<T> &mat)
{
  assert( (arr3d.dim2()==mat.num_rows())&(arr3d.dim3()==mat.num_cols()) );
  for (int i=0; i<mat.num_rows(); i++)
    for (int j=0; j<mat.num_cols(); j++)
      arr3d[position][i][j] = mat[i][j];
}

template <class T>
Matrix<T>  GetMatrixInArray3D(Array3D<T> &arr3d, int position)
{
  Matrix<double> mat(arr3d.dim2(), arr3d.dim3());
  for (int i=0; i<arr3d.dim2(); i++)
    for (int j=0; j<arr3d.dim3(); j++)
      mat[i][j] = arr3d[position][i][j];
  return mat;
}

template <class T>
static void   PutRowVectorInArray3D(Array3D<T> &arr3d, int posInDim1,
		int posInDim2, Vector<T> &vec)
{
  assert( arr3d.dim3() == vec.size() );
  for (int i=0; i<vec.size(); i++)
    arr3d[posInDim1][posInDim2][i] = vec[i];
}

template <class T>
Vector<T>  GetRowVectorInArray3D(Array3D<T> &arr3d, int posInDim1,
                int posInDim2)
{
  Vector<T>  vec(arr3d.dim3());
  for (int i=0; i<vec.size(); i++)
    vec[i] = arr3d[posInDim1][posInDim2][i];
  return  vec;
}


/* Put vector in  arr[ posInDim1, ..., posInDim3] */
template <class T>
static void   PutColVectorInArray3D(Array3D<T> &arr3d, int posInDim1,
                int posInDim3, Vector<T> &vec)
{
  assert( arr3d.dim2() == vec.size() );
  for (int i=0; i<vec.size(); i++)
    arr3d[posInDim1][i][posInDim3] = vec[i];
}

/* Get vector in  arr[ posInDim1, ..., posInDim3] */
template <class T>
Vector<T>  GetColVectorInArray3D(Array3D<T> &arr3d, int posInDim1,
                int posInDim3)
{
  Vector<T>  vec(arr3d.dim2());
  for (int i=0; i<vec.size(); i++)
    vec[i] = arr3d[posInDim1][i][posInDim3];
  return  vec;
}



/***********************************************************************
 ***      Math Functions
 ***********************************************************************/
/* Function: LogGamma
 * Usage: LogGamma(x)
 * -------------------
 * This function returns log Gamma function in double precision.
 */

double	LogGamma(double x)
{
	int	k;
	double	w, t, y, v;
	double	a[22] = {  9.967270908702825e-5, -1.9831672170162227e-4, -0.00117085315349625822, 0.00722012810948319552, 
			-0.0096221300936780297, -0.04219772092994235254, 0.16653861065243609743, -0.04200263501129018037, 
			-0.65587807152061930091, 0.57721566490153514421, 0.99999999999999999764, 
			4.67209725901142e-5, -6.812300803992063e-5, -0.00132531159076610073, 0.0073352117810720277, 
			-0.00968095666383935949, -0.0421764281187354028, 0.16653313644244428256, -0.04200165481709274859, 
			-0.65587818792782740945, 0.57721567315209190522, 0.99999999973565236061 };
	double b[98] = { -4.587497028e-11, 1.902363396e-10, 8.6377323367e-10, 1.15513678861e-8, -2.556403058605e-8, -1.5236723372486e-7, 
        		-3.1680510638574e-6, 1.22903704923381e-6, 2.334372474572637e-5, 0.00111544038088797696, 
			0.00344717051723468982, 0.03198287045148788384, -0.32705333652955399526, 0.40120442440953927615, 
			-5.184290387e-11, -8.3355121068e-10, -2.56167239813e-9, 1.455875381397e-8, 1.3512178394703e-7, 2.9898826810905e-7, 
			-3.58107254612779e-6, -2.445260816156224e-5, -4.417127762011821e-5, 0.00112859455189416567, 
			0.00804694454346728197, 0.04919775747126691372, -0.24818372840948854178, 0.11071780856646862561, 
			3.0279161576e-10, 1.60742167357e-9, -4.05596009522e-9, -5.089259920266e-8, -2.029496209743e-8, 1.35130272477793e-6, 
			3.91430041115376e-6, -2.871505678061895e-5, -2.3052137536922035e-4, 4.5534656385400747e-4, 
			0.01153444585593040046, 0.07924014651650476036, -0.12152192626936502982, -0.07916438300260539592, 
			-5.091914958e-10, -1.15274986907e-9, 1.237873512188e-8, 2.937383549209e-8, -3.0621450667958e-7, -7.7409414949954e-7, 
			8.16753874325579e-6, 2.412433382517375e-5, -2.60612176060637e-4, -9.1000087658659231e-4, 
			0.01068093850598380797, 0.11395654404408482305, 0.07209569059984075595, -0.10971041451764266684, 
			4.0119897187e-10, -1.3224526679e-10, -1.002723190355e-8, 2.569249716518e-8, 2.0336011868466e-7, -1.1809768272606e-6, 
			-3.00660303810663e-6, 4.402212897757763e-5, -1.462405876235375e-5, -0.0016487379559600128, 
			0.00513927520866443706, 0.13843580753590579416, 0.32730190978254056722, 0.08588339725978624973, 
			-1.5413428348e-10, 6.4905779353e-10, 1.60702811151e-9, -2.655645793815e-8, 7.619544277956e-8, 4.7604380765353e-7, 
			-4.90748870866195e-6, 8.21513040821212e-6, 1.4804944070262948e-4, -0.00122152255762163238, 
			-8.7425289205498532e-4, 0.1443870369965796831, 0.61315889733595543766, 0.55513708159976477557, 
			1.049740243e-11, -2.5832017855e-10, 1.39591845075e-9, -2.1177278325e-10, -5.082950464905e-8, 3.7801785193343e-7, 
			-7.3982266659145e-7, -1.088918441519888e-5, 1.2491810452478905e-4, -4.9171790705139895e-4, 
			-0.0042570708944826646, 0.13595080378472757216, 0.89518356003149514744, 1.31073912535196238583 };
    double c[65] = {
        1.16333640008e-8, -8.33156123568e-8, 3.832869977018e-7, -1.5814047847688e-6,  6.50106723241e-6, -2.74514060128677e-5, 
        1.209015360925566e-4, -5.666333178228163e-4, 0.0029294103665559733, -0.0180340086069185819, 0.1651788780501166204, 
	1.1031566406452431944, 1.2009736023470742248, 1.3842760642e-9, -6.9417501176e-9, 3.42976459827e-8, -1.785317236779e-7, 
        9.525947257118e-7, -5.2483007560905e-6, 3.02364659535708e-5, -1.858396115473822e-4, 
        0.0012634378559425382, -0.0102594702201954322, 0.1243625515195050218, 1.3888709263595291174, 2.4537365708424422209, 
        1.298977078e-10, -8.02957489e-10, 4.945484615e-9, -3.17563534834e-8, 2.092136698089e-7, -1.4252023958462e-6, 
        1.01652510114008e-5, -7.74550502862323e-5, 6.537746948291078e-4, -0.006601491253552183, 0.0996711934948138193, 
	1.6110931485817511402, 3.9578139676187162939, 1.83995642e-11, -1.353537034e-10, 9.984676809e-10, -7.6346363974e-9, 
        5.99311464148e-8, -4.868554120177e-7, 4.1441957716669e-6, -3.77160856623282e-5, 3.805693126824884e-4, -0.0045979851178130194, 
        0.0831422678749791178, 1.7929113303999329439, 5.6625620598571415285, 3.4858778e-12, -2.97587783e-11, 
        2.557677575e-10, -2.2705728282e-9, 2.0702499245e-8, -1.954426390917e-7, 1.9343161886722e-6, -2.0479024910257e-5, 
        2.405181940241215e-4, -0.0033842087561074799, 0.0713079483483518997, 1.9467574842460867884, 7.5343642367587329552
    };
    double d[7] = {
        -0.00163312359200500807, 8.3644533703385956e-4, -5.9518947575728181e-4, 7.9365057505415415e-4, 
        -0.00277777777735463043, 0.08333333333333309869, 0.91893853320467274178 };

    w = x;
    if (x < 0) {
        w = 1 - x;
    }
    if (w < 0.5) {
        k = w < 0.25 ? 0 : 11;
        y = ((((((((((a[k] * w + a[k + 1]) * w + a[k + 2]) * w + a[k + 3]) * w + a[k + 4]) * w + 
            a[k + 5]) * w + a[k + 6]) * w + a[k + 7]) * w + a[k + 8]) * w + a[k + 9]) * w + a[k + 10]) * w;
        y = -log(y);
    } else if (w < 3.5) {
        t = w - 4.5 / (w + 0.5);
        k = (int) (t + 4);
        t -= k - 3.5;
        k *= 14;
        y = ((((((((((((b[k] * t + b[k + 1]) * t + b[k + 2]) * t + b[k + 3]) * t + b[k + 4]) * t + 
            b[k + 5]) * t + b[k + 6]) * t + b[k + 7]) * t + b[k + 8]) * t + b[k + 9]) * t + b[k + 10]) * t + 
            b[k + 11]) * t + b[k + 12]) * t + b[k + 13];
    } else if (w < 8) {
        k = ((int) w) - 3;
        t = w - (k + 3.5);
        k *= 13;
        y = (((((((((((c[k] * t + c[k + 1]) * t + c[k + 2]) * t + c[k + 3]) * t + c[k + 4]) * t + 
            c[k + 5]) * t + c[k + 6]) * t + c[k + 7]) * t + c[k + 8]) * t + c[k + 9]) * t + c[k + 10]) * t + 
            c[k + 11]) * t + c[k + 12];
    } else {
        v = 1 / w;
        t = v * v;
        y = (((((d[0] * t + d[1]) * t + d[2]) * t + d[3]) * t + d[4]) * t + d[5]) * v + d[6];
        y += (w - 0.5) * log(w) - w;
    }
    if (x < 0) {
        y = log(3.141592653589793238 / sin(x * 3.141592653589793238)) - y;
    }
    return y;
}

double  Density_StudentT( double x, double v, double mu, double sigma2)
{
  return ( pow( 2.7183, LogGamma(v/2.0 + 0.5) - LogGamma(0.5*v) )
       / sqrt( sigma2*v*3.1416 ) * pow( 1 + (x-mu)*(x-mu)/sigma2/v , 
       -0.5*v-0.5 ) );
}

double   Multi_NormalDen(Vector<double>  x, Vector<double>  mu,
			 Matrix<double>  cov)
{
  assert( (mu.dim()==cov.num_rows()) && (mu.dim()==cov.num_cols() ));
  int    dim=mu.dim();
  Vector<double>   diff = x - mu;
  double logden = - 0.5*log(2*3.1416) - 0.5*log(det(cov))
                  - 0.5*mul_vecmatvec(diff, cov);
  return pow(2.7183, logden);
}



template <class T>
T		Sum(const Vector<T>   &vec)
{
  T   sum=vec[0]-vec[0];
  for (int i=0; i<vec.size(); i++)	sum += vec[i];
  return  sum;
}

template <class T>
T         Max(const Vector<T>   &vec)
{
  T      tmp=vec[0];
  for (int i=0; i<vec.dim(); i++) tmp=(tmp<=vec[i])? vec[i]:tmp;
  return       tmp;
}

template <class T>
T         Min(T  a, T  b, T c)
{
  T       tmp = (a<b)? a:b;
  return  ((tmp<c)?tmp:c);
}



/* This function normalize the vector */
static void  NormalizeWeights(Vector<double> &weights, int length)
{
  assert( weights.size() >= length );
  double total=0.0;
  for (int i=0; i<length; i++) total += weights[i];
  for (int i=0; i<length; i++) weights[i] = weights[i]/total;
}

/* This function returns the index of the mimimum element in the first 'numEle' 
elements. */              
static int  IndexofMinValue( Vector<double> &vec, int  numEle)
{
  assert( vec.size() >= numEle );
  int     index=0;
  double  minV = vec[0];
  for ( int i = 0; i < numEle ; i++)
    if ( minV > vec[i] ) {  minV = vec[i];	index = i;  }
  return  index;
} 

/* Notice that startPos and endPos here are actual positions. */
static int  CountSymbol( const Vector<int> &vec, int startPos,
                int endPos, int symbol)
{
  int   count=0;
  for (int i=startPos; i<=endPos; i++)
    count += ( ( vec[i-1] == symbol ) ? 1:0 );
  return  count;
}

/* Compute the empirical CV */
double	GetCV2(Vector<double> vec)
{
  double   var=0.0, mean = Sum(vec)/vec.size();
  for (int i=0; i<vec.size(); i++)  var += (vec[i]-mean)*(vec[i]-mean);
  return   (var/(vec.size()-1)/mean/mean);
}





#endif

