static const char *matrixMultiply = 
"float4 main(in float2 position : TEX0,                             \n"
"            uniform int  startOps,                                 \n"
"            uniform int  stopOps,                                  \n"
"            uniform samplerRECT  accum,                            \n"
"            uniform samplerRECT  texLHS,                           \n"
"			 uniform samplerRECT  texRHS) : COLOR                   \n"
"{																	\n"
"	 float4 output = texRECT(accum, position); 						\n"
"	 for (int i = startOps; i < stopOps; i++) {						\n"
"		float4 a  = texRECT(texLHS, float2(i,position.y));          \n"
"		float4 b  = texRECT(texRHS, float2(position.x,i));          \n"
"		output += a.xxzz*b.xyxy;					                \n"
"		output += a.yyww*b.zwzw;									\n"
"    }																\n"
"    return output;                                                 \n"
"}                                                                  \n";

static const char *testInputs = 
"float4 main(in float2 position : TEX0,                             \n"
"            uniform int  startOps,                                 \n"
"            uniform int  stopOps,                                  \n"
"            uniform samplerRECT  accum,                            \n"
"            uniform samplerRECT  texLHS,                           \n"
"			 uniform samplerRECT  texRHS) : COLOR                   \n"
"{																	\n"
"	 float4 a  = texRECT(texLHS, position);						    \n"
"	 float4 b  = texRECT(texRHS, position);							\n"
"	 float4 c  = float4(startOps, stopOps,position);				\n"
"    return a;										                \n"
"}                                                                  \n";

/*
"		output.x += a.x*b.x;					            \n"
"		output.x += a.y*b.z;					            \n"
"		output.y += a.x*b.y;					            \n"
"		output.y += a.y*b.w;					            \n"
"		output.z += a.z*b.x;					            \n"
"		output.z += a.w*b.z;					            \n"
"		output.w += a.z*b.y;					            \n"
"		output.w += a.w*b.w;					            \n"
*/

static const char *addConstant = 
"float4 main(in float2 coords : TEX0,                               \n"
"            uniform samplerRECT  texture,                          \n"
"			 uniform float value) : COLOR                           \n"
"{                                                                  \n"       
"    float4 c  = texRECT(texture, coords);                          \n"
"    return c + value;                                              \n"
"}                                                                  \n";

static const char *multiplyConstant = 
"float4 main(in float2 coords : TEX0,                               \n"
"            uniform samplerRECT  texture,                          \n"
"			 uniform float value) : COLOR                           \n"
"{                                                                  \n"
"    float4 c  = texRECT(texture, coords);                          \n"
"    return c * value;                                              \n"
"}                                                                  \n";

static const char *addTextures = 
"float4 main(in float2 coords : TEX0,                               \n"
"            uniform samplerRECT  texLHS,                           \n"
"			 uniform samplerRECT  texRHS) : COLOR                   \n"
"{                                                                  \n"
"    float4 l  = texRECT(texLHS, coords);                           \n"
"    float4 r  = texRECT(texRHS, coords);                           \n"
"    return l + r;                                                  \n"
"}                                                                  \n";

static const char *copyTexture = 
"float4 main(in float2 coords : TEX0,                               \n"
"            uniform samplerRECT  texture) : COLOR                  \n"
"{                                                                  \n"       
"    return texRECT(texture, coords);                               \n"
"}                                                                  \n";

//doesn't take fringes into accout... this shader is fine if the data is even dimensioned
static const char *assignConst = 
"float4 main(uniform float v) : COLOR                               \n"
"{                                                                  \n"
"    return v;                                                      \n"
"}                                                                  \n";

static const char *fixRightFringe = 
"float4 main(in float2 coords : TEX0,                               \n"
"			 uniform float value,                                   \n"
"            uniform samplerRECT  texture) : COLOR                  \n"
"{                                                                  \n"
"    float4 v  = texRECT(texture, coords);                          \n"
"	 if(coords.x > value){                                          \n"
"        return float4(v.x, 0, v.z, 0);                             \n"
"    } else {                                                       \n"
"        return v;                                                  \n"
"    }                                                              \n"
"}                                                                  \n";

static const char *fixBottomFringe = 
"float4 main(in float2 coords : TEX0,                               \n"
"			 uniform float value,                                   \n"
"            uniform samplerRECT  texture) : COLOR                  \n"
"{                                                                  \n"
"    float4 v  = texRECT(texture, coords);                          \n"
"	 if(coords.y > value){                                          \n"
"        return float4(v.x, v.y, 0, 0);                             \n"
"    } else {                                                       \n"
"        return v;                                                  \n"
"    }                                                              \n"
"}                                                                  \n";