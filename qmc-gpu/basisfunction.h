#include <string>
#include "matrix.h"

typedef GLfloat finalType;

static const char * shaderHead =
"float4 main(in float2 coords : TEX0,                                   \n"
"            uniform samplerRECT  epos) : COLOR                         \n"
"{                                                                      \n"
"   float3 r  = texRECT(epos, coords).xyz;                              \n"
"   float3 klm = float3(%i,%i,%i);                                      \n"
"   float3 center = float3(%E,%E,%E);                                   \n"
"   float ngauss = %i;                                                  \n"
"   float4 output = float4(0,0,0,0);                                    \n"
"   float r_sq = 0;                                                     \n"
"   float xyz_term = 0;                                                 \n"
"   float coeffs[%i][2] = {                                             \n";

static const char * shaderCoeff =
"       {%E, %E},                                                       \n";

static const char * shaderTail =
"   };                                                                  \n"
"   r = r - center;                                                     \n"
"   r_sq = dot(r,r);                                                    \n"
"   xyz_term = pow(r.x,klm.x)*pow(r.y,klm.y)*pow(r.z,klm.z);            \n" 
"   for(int j=0; j<ngauss; j++){                                        \n"
"       output.x += coeffs[j][1]*exp(-coeffs[j][0]*r_sq);               \n"
"   }                                                                   \n"
"   output.x *= xyz_term;                                               \n"
"                                                                       \n"
"   float exp_xyz_term = 0;                                             \n"
"   float temp = 0;                                                     \n"
"   for(int j=0; j<ngauss; j++){                                        \n"
"       exp_xyz_term = coeffs[j][1]*exp(-coeffs[j][0]*r_sq)*xyz_term;   \n"
"       temp = -2.0*coeffs[j][0];                                       \n"
"                                                                       \n"
"       output.y += (klm.x/r.x + temp*r.x)*exp_xyz_term;                \n"
"       output.z += (klm.y/r.y + temp*r.y)*exp_xyz_term;                \n"
"       output.w += (klm.z/r.z + temp*r.z)*exp_xyz_term;                \n"
"   }                                                                   \n"
"   return output;                                                      \n"
"}                                                                      \n";

/**
the fundamental ideology behind the design of this class is this object is set
for a given set of parameters (the coefficients etc). after one of these objects
has been instantiated, it is capable of calculating basis functions for
any set of electronic positions. e.g. it is a streaming function that operates:
R >> calculateBasisFunction >> psi, gradient, laplacian

perhaps this will need to be modified to:
R, Rc >> calculateBasisFunction >> psi, gradient, laplacian

where R and Rc are vectors of electronic and nuclear coordinates respectively.

for this to be moved to QMcBeaver, i need to figure out for a given number
of these BasisFunction objects, how many walkers can i run at a time for a given
amount of memory? (ram and/or graphics memory)

if all the data can be kept "online" then other shaders (if appropriate) can operate on the data.

lastly, if data needs to be removed from the graphics card, then there needs to be implemented
a more global way of keeping track of the texture ids used.
*/
class BasisFunction : public Matrix {
public:	
	/*basically just the info for calculating a basis function. my idea is that 
	eventually, most/all of these parameters would be "hard coded" into a char *
	just like what is done in generateBasisFunction.*/
	BasisFunction(int _numWalkers, Array2D<double> _Coeffs, int _K, int _L, int _M, 
                  int _N_Gauss, double _Xc, double _Yc, double _Zc, string _Type){
		dim = sqrt((double)_numWalkers);//require a perfect square for now
		assert( dim*dim == _numWalkers );
		setData(_numWalkers,_Coeffs,_K,_L,_M,_N_Gauss,_Xc, _Yc, _Zc, _Type);
	}

	bool calculateBasisFunctions(Array2D<double> &R, bool useGPU){
		assert( R.dim1() == numWalkers );
		calculated = true;
		if(useGPU){
			calculateWithGPU(R);
		} else {
			calculateWithCPU(R);
		}
		return calculated;
	}

	void getPsi(const int whichBF, finalType * psi){
		if(!calculated){
			printf("calculate it first!\n");
			return;
		}
		(*psi) = pixelData[whichBF];
	}

	void getGradientPsi(const int whichBF, finalType * gradPsiX, finalType * gradPsiY, finalType * gradPsiZ){
		if(!calculated){
			printf("calculate it first!\n");
			return;
		}
		(*gradPsiX) = pixelData[whichBF +1];
		(*gradPsiY) = pixelData[whichBF +2];
		(*gradPsiZ) = pixelData[whichBF +3];
	}

private:
	void calculateWithGPU(Array2D<double> &R){		
		int i, j, index;
		for(i=0; i<dim; i++){
            for(j=0; j<dim; j++){
				index = mapping(i, j, dim, dim);
				pixelData[index   ] = R(i*dim+j,0);
                pixelData[index +1] = R(i*dim+j,1);
                pixelData[index +2] = R(i*dim+j,2);
                pixelData[index +3] = 0;
            }
        }		

		GLuint texId[1];
		glGenTextures(1, texId);	
		glActiveTextureARB(GL_TEXTURE0_ARB);
		glBindTexture(GL_TEXTURE_RECTANGLE_NV, texId[0]);
		glTexImage2D(GL_TEXTURE_RECTANGLE_NV, 0, GL_FLOAT_RGBA32_NV, 
					 dim, dim, 0, GL_RGBA, GL_FLOAT, pixelData);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		
		textureData->BeginCapture();
		
		cgGLEnableProfile(g_cgProfile);	
		cgGLBindProgram(fragProg);
		cgGLSetTextureParameter(tex, texId[0]);
		cgGLEnableTextureParameter(tex);

		glClearColor(0.0f, 0.05f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		int maxs = textureData->GetMaxS();
		int maxt = textureData->GetMaxT(); 
		glBegin(GL_QUADS);
			glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
			glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
			glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
			glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
		glEnd();

		textureData->EndCapture();

		glFinish();
		glFlush();
		glDeleteTextures(1, texId);
		getError("Error in calculateWithGPUText function");	
		unloadMatrix(false);
	}

	//a way needs to be created for this method to store data in double precision
	void calculateWithCPU(Array2D<double> &R){
		//data = Array2D<finalType>(2*dim,2*dim);
		//finalType** collection = data.array();

		double xyz_term, r_sq, radialFunction=0;
		double x,y,z, gx=0, gy=0, gz=0;
		int index;
		for(int i=0; i<numWalkers; i++){
			x = R(i,0) - xc;
			y = R(i,1) - yc;
			z = R(i,2) - zc;
			r_sq = x*x + y*y + z*z;
			xyz_term = pow(x,k)*pow(y,l)*pow(z,m);

			radialFunction = 0.0;
			gx = 0.0;
			gy = 0.0;
			gz = 0.0;
			for(int j=0; j<N_Gauss; j++){
				//first calculate psi		
				radialFunction += Coeffs(j,1)*exp(-Coeffs(j,0)*r_sq);

				//2nd, calculate the gradient
				double exp_xyz_term = Coeffs(j,1)*exp(-Coeffs(j,0)*r_sq)*xyz_term;
				double temp = -2.0*Coeffs(j,0);

				gx += (k/x + temp*x)*exp_xyz_term;
				gy += (l/y + temp*y)*exp_xyz_term;
				gz += (m/z + temp*z)*exp_xyz_term;
			}

			index = mapping((int)(i/dim),(i%dim),dim,dim);
			pixelData[index  ] = xyz_term * radialFunction;
			pixelData[index+1] = gx;
			pixelData[index+2] = gy;
			pixelData[index+3] = gz;
		}
	}

	void setData(int _numWalkers, Array2D<double> _Coeffs, int _K, int _L, int _M, 
		int _N_Gauss, double _Xc, double _Yc, double _Zc, string _Type){
		Coeffs = _Coeffs;
		k = _K; l = _L; m = _M;
		N_Gauss = _N_Gauss;
		Type = _Type;
		numWalkers = _numWalkers;
		calculated = false;
		xc = _Xc; yc = _Yc; zc = _Zc;
		generateShader(_Coeffs,_K,_L,_M,_Xc,_Yc,_Zc,_N_Gauss);
		fragProg = cgCreateProgram(g_cgContext, CG_SOURCE,
                   shader, g_cgProfile,
                   "main", NULL);		
		cgGLLoadProgram(fragProg);		
		tex = cgGetNamedParameter(fragProg, "epos");
		textureData = new RenderTexture(textureMode);
		textureData->Initialize(dim,dim,true,false);
		rows = dim*2;
		columns = dim*2;
		pixelData = new GLfloat[dim*dim*4];
	}

	void generateShader(Array2D<double> & coeff,
						const int k, const int l, const int m,
						const double xc, const double yc, const double zc,
						const int ngauss){
		char temp[2048];
		shader = (char*)malloc(100*1024);
		sprintf(shader, shaderHead,k,l,m,xc,yc,zc,ngauss,ngauss);
		for(int i=0; i<ngauss; i++){
			sprintf(temp,shaderCoeff,coeff(i,0),coeff(i,1));
			//i think this next line assumes that temp always has the same # of chars
			strcat(shader, temp);
		}
		strcat(shader,shaderTail);
		//cout << shader << endl;
	}

	int numWalkers, dim;
	double xc, yc, zc;
	CGprogram fragProg;
	CGparameter tex;
	char * shader;
	bool calculated;
	
	/**
    Array containing the parameters for the basis functions where
    Coeffs[Gaussian #][0=exp,1=contract]
    */
	Array2D <double> Coeffs;

	/**
    Array containing the k,l,m parameters which indicate the 
    "angular momentum state" of the basis function 
    (\f$bf=x^{k}y^{l}z^{m}*RadialFunction(r)\f$) where xyz[bf #][0=k,1=l,2=m].
    For example, a "px" orbital would have \f$(k,l,m)=(1,0,0)\f$.
    */
	int k,l,m;                       
                           
	/**
    Array containing the number of gaussians that need to be contracted
    for the radial portion of the basis function 
    (\f$bf=x^{k}y^{l}z^{m}*RadialFunction(r)\f$) where N_Gauss[bf #].
    */
	int N_Gauss;                                
                              
	/**
    Array containing the type of the basis function where Type[bf #].
    The type is a string representation of the "angular momentum state."
    For example, "px", "dxy", and "fxxx" are all types of basis functions.
    */
	string Type; 
};