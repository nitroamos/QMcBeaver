#include <string>
#include "matrix.h"

typedef GLfloat finalType;

static const char *generateBasisFunction = 
"float4 main(in float3 r : TEX7) : COLOR                             \n"
"{																	\n"
"	float3 klm = float3(0,0,0);		    							\n"
"	float3 center = float3(0,0,0);         							\n"
"	float ngauss = 8;												\n"
"	float4 output = float4(0,0,0,0);         						\n"
"	float r_sq = 0;													\n"
"	float xyz_term = 0;												\n"
"	float coeffs[8][2] = {											\n"
"		{6.66500000E+03,	0.363584299905},                        \n"
"		{1.00000000E+03,	0.674985792448},						\n"
"		{2.28000000E+02,	1.1316199392},							\n"
"		{6.47100000E+01,	1.65300920982},							\n"
"		{2.10600000E+01,	1.92382064638},							\n"
"		{7.49500000E+00,	1.44727831645},							\n"
"		{2.79700000E+00,	0.439163129646},						\n"
"		{5.21500000E-01,	0.00664580366553}};						\n"
"	r = r - center;													\n"
"	r_sq = dot(r,r);												\n"
"	xyz_term = pow(r.x,klm.x)*pow(r.y,klm.y)*pow(r.z,klm.z);		\n"	
"	for(int j=0; j<ngauss; j++){									\n"
"		output.x += coeffs[j][1]*exp(-coeffs[j][0]*r_sq);			\n"
"	}																\n"
"	output.x *= xyz_term;											\n"
"																	\n"
"	float exp_xyz_term = 0;											\n"
"	float temp = 0;													\n"
"	for(int j=0; j<ngauss; j++){									\n"
"		exp_xyz_term = coeffs[j][1]*exp(-coeffs[j][0]*r_sq)*xyz_term; \n"
"		temp = -2.0*coeffs[j][0];									\n"
"																	\n"
"		output.y += (klm.x/r.x + temp*r.x)*exp_xyz_term;			\n"
"		output.z += (klm.y/r.y + temp*r.y)*exp_xyz_term;			\n"
"		output.w += (klm.z/r.z + temp*r.z)*exp_xyz_term;			\n"
"	}																\n"
"   return output;	                                                \n"
"}                                                                  \n";

static const char *generateBasisFunctionTexture = 
"float4 main(in float2 coords : TEX0,                               \n"
"			 uniform samplerRECT  epos) : COLOR                     \n"
"{																	\n"
"   float3 r  = texRECT(epos, coords).xyz;	                        \n"
"	float3 klm = float3(0,0,0);		    							\n"
"	float3 center = float3(0,0,0);         							\n"
"	float ngauss = 8;												\n"
"	float4 output = float4(0,0,0,0);         						\n"
"	float r_sq = 0;													\n"
"	float xyz_term = 0;												\n"
"	float coeffs[8][2] = {											\n"
"		{6.66500000E+03,	0.363584299905},                        \n"
"		{1.00000000E+03,	0.674985792448},						\n"
"		{2.28000000E+02,	1.1316199392},							\n"
"		{6.47100000E+01,	1.65300920982},							\n"
"		{2.10600000E+01,	1.92382064638},							\n"
"		{7.49500000E+00,	1.44727831645},							\n"
"		{2.79700000E+00,	0.439163129646},						\n"
"		{5.21500000E-01,	0.00664580366553}};						\n"
"	r = r - center;													\n"
"	r_sq = dot(r,r);												\n"
"	xyz_term = pow(r.x,klm.x)*pow(r.y,klm.y)*pow(r.z,klm.z);		\n"	
"	for(int j=0; j<ngauss; j++){									\n"
"		output.x += coeffs[j][1]*exp(-coeffs[j][0]*r_sq);			\n"
"	}																\n"
"	output.x *= xyz_term;											\n"
"																	\n"
"	float exp_xyz_term = 0;											\n"
"	float temp = 0;													\n"
"	for(int j=0; j<ngauss; j++){									\n"
"		exp_xyz_term = coeffs[j][1]*exp(-coeffs[j][0]*r_sq)*xyz_term; \n"
"		temp = -2.0*coeffs[j][0];									\n"
"																	\n"
"		output.y += (klm.x/r.x + temp*r.x)*exp_xyz_term;			\n"
"		output.z += (klm.y/r.y + temp*r.y)*exp_xyz_term;			\n"
"		output.w += (klm.z/r.z + temp*r.z)*exp_xyz_term;			\n"
"	}																\n"
"   return output;	                                                \n"
"}                                                                  \n";

static const char *inputCheck = 
"float4 main(in float3 r : TEX7) : COLOR							\n"
"{																	\n"
"   return float4(r,0);										        \n"
"}                                                                  \n";

static const char *simpleVertex = 
"void main(	float4 vpos				: POSITION,						\n"
"			float3 epos				: NORMAL,						\n"
"	        out float4 HPosition    : POSITION,						\n"
"			uniform float4x4 ModelViewProj,							\n"
"	        out float3 oe      : TEXCOORD7)		    				\n"
"{																	\n"
"	HPosition = mul(ModelViewProj, vpos);							\n"
//"	HPosition = float4(vpos,0,0);									\n"
//"   oe = float4(vpos.xy,epos.x,0);					                \n"
"   oe = epos;										                \n"
"}                                                                  \n";


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
			//calculateWithGPUVertex(R);
			calculateWithGPUText(R);
		} else {
			calculateWithCPU(R);
		}
		return calculated;
	}

	void getPsi(const int whichBF, finalType * psi){
		if(!calculated || !data.array()){
			printf("calculate it first!\n");
			return;
		}
		(*psi) = data( 2*(int)(whichBF/dim), 2*(whichBF%dim) );
	}

	void getGradientPsi(const int whichBF, finalType * gradPsiX, finalType * gradPsiY, finalType * gradPsiZ){
		if(!calculated){
			printf("calculate it first!\n");
			return;
		}
		if(!data.array()) return;
		(*gradPsiX) = data( 2*(int)(whichBF/dim) +1, 2*(whichBF%dim)   );
		(*gradPsiY) = data( 2*(int)(whichBF/dim)   , 2*(whichBF%dim) +1);
		(*gradPsiZ) = data( 2*(int)(whichBF/dim) +1, 2*(whichBF%dim) +1);
	}

private:
	void calculateWithGPUText(Array2D<double> &R){
		Stopwatch sw = Stopwatch();
		GLfloat* texels = new GLfloat[dim*dim*4];
		GLuint    texId[1];

		int i, j, index;
		for(i=0; i<dim; i++){
            for(j=0; j<dim; j++){
				index = mapping(i, j, dim, dim);
				texels[index   ] = R(i*dim+j,0);
                texels[index +1] = R(i*dim+j,1);
                texels[index +2] = R(i*dim+j,2);
                texels[index +3] = 0;
            }
        }
		

		glGenTextures(1, texId);	

		glActiveTextureARB(GL_TEXTURE0_ARB);
		glBindTexture(GL_TEXTURE_RECTANGLE_NV, texId[0]);
		sw.reset();
		sw.start();
		
		glTexImage2D(GL_TEXTURE_RECTANGLE_NV, 0, GL_FLOAT_RGBA32_NV, 
							dim, dim, 0, GL_RGBA, GL_FLOAT, texels);
		sw.stop();
		cout << "total uploading time " << sw.timeMS() << endl;
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_S, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_WRAP_T, GL_CLAMP);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameterf(GL_TEXTURE_RECTANGLE_NV, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

		CGprogram fragProg;
		CGparameter tex;
		if(true){
			fragProg = cgCreateProgram(g_cgContext, CG_SOURCE,
                   generateBasisFunctionTexture, g_cgProfile,
                   "main", NULL);
		} else {
			fragProg = cgCreateProgram(g_cgContext, CG_SOURCE,
                   inputCheck, g_cgProfile,
                   "main", NULL);
		}
		
		if(fragProg != NULL){
			cgGLLoadProgram(fragProg);
			tex = cgGetNamedParameter(fragProg, "epos");
		}

		RenderTexture * result = new RenderTexture(textureMode);
		result->Initialize(dim,dim,true,false);
		
		result->BeginCapture();
		
		cgGLEnableProfile(g_cgProfile);	
		cgGLBindProgram(fragProg);
		
		cgGLSetTextureParameter(tex, texId[0]);
		cgGLEnableTextureParameter(tex);


		glClearColor(0.0f, 0.05f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		
		int maxs = result->GetMaxS();
		int maxt = result->GetMaxT(); 
		glBegin(GL_QUADS);
			glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
			glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
			glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
			glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
		glEnd();

		result->EndCapture();
		
		glFinish();
		glFlush();

		getError("Error in calculateWithGPUText function");
		if(textureData) delete textureData;
		textureData = result;
		rows = dim*2;
		columns = dim*2;
		sw.reset();
		sw.start();
		unloadMatrix(false);
		sw.stop();
		cout << "total unloading time " << sw.timeMS() << endl;
	}

	void calculateWithGPUVertex(Array2D<double> &R){
		int texDim = dim;
		GLfloat * vertPos = new GLfloat[2*numWalkers];
		GLfloat * elecPos=  new GLfloat[3*numWalkers];

		double temp;
		for(int i=0; i<numWalkers; i++){
			temp = (int)(i/dim);
			vertPos[2*i+1] = -1 + 2*temp/(dim-1);
			temp = i%dim;
			vertPos[2*i  ] = -1 + 2*temp/(dim-1);
			if(true){
				double shift = 0.9999;//what a screwy solution
				vertPos[2*i  ] *= shift;
				vertPos[2*i+1] *= shift;
			}
			elecPos[3*i  ] = R(i,0);
			elecPos[3*i+1] = R(i,1);
			elecPos[3*i+2] = R(i,2);
			//printf("(%4.2g, %4.2g) is (%8.5g, %8.5g, %8.5g)\n",vertPos[2*i],vertPos[2*i+1],elecPos[3*i],elecPos[3*i+1],elecPos[3*i+2]);	
		}

		CGprogram vertProg, fragProg;
		CGprofile vertProfile = cgGLGetLatestProfile(CG_GL_VERTEX);

		vertProg = cgCreateProgram(g_cgContext, CG_SOURCE,
	               simpleVertex, vertProfile,
		           "main", NULL);

		if(vertProg != NULL) cgGLLoadProgram(vertProg);

		if(true){
			fragProg = cgCreateProgram(g_cgContext, CG_SOURCE,
                   generateBasisFunction, g_cgProfile,
                   "main", NULL);
		} else {
			fragProg = cgCreateProgram(g_cgContext, CG_SOURCE,
                   inputCheck, g_cgProfile,
                   "main", NULL);
		}

        if(fragProg != NULL) cgGLLoadProgram(fragProg);

		RenderTexture * result = new RenderTexture(textureMode);
		result->Initialize(texDim,texDim,true,false);
		
		result->BeginCapture();
		
		cgGLEnableProfile(g_cgProfile);	
		cgGLBindProgram(fragProg);
			
		glClearColor(0.0f, 0.05f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		if(false){
			int maxs = result->GetMaxS();
			int maxt = result->GetMaxT(); 
			glBegin(GL_QUADS);
				glTexCoord2f(0.0,  0.0);   glVertex3f(-1.0, -1.0, 0.0);
				glTexCoord2f(0.0,  maxt);  glVertex3f(-1.0, 1.0, 0.0);
				glTexCoord2f(maxs, maxt);  glVertex3f(1.0, 1.0, 0.0);
				glTexCoord2f(maxs, 0.0);   glVertex3f(1.0, -1.0, 0.0);
			glEnd();
		} else {
			cgGLEnableProfile(vertProfile);
			cgGLBindProgram(vertProg);
			cgGLSetStateMatrixParameter(cgGetNamedParameter(vertProg, "ModelViewProj"),
			CG_GL_MODELVIEW_PROJECTION_MATRIX,
            CG_GL_MATRIX_IDENTITY);
			CGparameter vpos = cgGetNamedParameter(vertProg, "vpos");
			cgGLEnableClientState(vpos);
			cgGLSetParameterPointer(vpos, 2, GL_FLOAT, 0, vertPos);		
			CGparameter epos = cgGetNamedParameter(vertProg, "epos");
			cgGLEnableClientState(epos);
			cgGLSetParameterPointer(epos, 3, GL_FLOAT, 0, elecPos);

			glDrawArrays(GL_POINTS,0,numWalkers);

			cgGLDisableClientState(vpos);
			cgGLDisableClientState(epos);
		}

		result->EndCapture();
		
		glFinish();
		glFlush();
		
		cgGLDisableProfile(g_cgProfile);
		cgGLDisableProfile(vertProfile);

		getError("Error in calculateWithGPUVertex function");
		if(textureData) delete textureData;
		textureData = result;
		rows = texDim*2;
		columns = texDim*2;
		unloadMatrix(false);
	}

	void calculateWithCPU(Array2D<double> &R){
		data = Array2D<finalType>(2*dim,2*dim);
		finalType** collection = data.array();
		
		double xyz_term, r_sq, radialFunction=0;
		double **coeffs;
		double x,y,z, gx=0, gy=0, gz=0;
		for(int i=0; i<numWalkers; i++){
			x = R(i,0) - xc;
			y = R(i,1) - yc;
			z = R(i,2) - zc;
			r_sq = x*x + y*y + z*z;
			xyz_term = pow(x,k)*pow(y,l)*pow(z,m);

			//first calculate psi		
			radialFunction = 0.0;
			coeffs = Coeffs.array();

			for(int j=0; j<N_Gauss; j++){
				//cout << "coeffs are " << coeffs[j][0] << " and " << coeffs[j][1] << endl;
				radialFunction += coeffs[j][1]*exp(-coeffs[j][0]*r_sq);
			}

			collection[ 2*(int)(i/dim) ][ 2*(i%dim) ] = xyz_term * radialFunction;

			//2nd, calculate the gradient
			gx = 0.0;
			gy = 0.0;
			gz = 0.0;

			for(int j=0; j<N_Gauss; j++){
				double exp_xyz_term = coeffs[j][1]*exp(-coeffs[j][0]*r_sq)*xyz_term;
				double temp = -2.0*coeffs[j][0];

				gx += (k/x + temp*x)*exp_xyz_term;
				gy += (l/y + temp*y)*exp_xyz_term;
				gz += (m/z + temp*z)*exp_xyz_term;
			}

			collection[ 2*(int)(i/dim) +1][ 2*(i%dim)   ] = gx;
			collection[ 2*(int)(i/dim)   ][ 2*(i%dim) +1] = gy;
			collection[ 2*(int)(i/dim) +1][ 2*(i%dim) +1] = gz;
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
	}

	int numWalkers, dim;
	double xc, yc, zc;

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