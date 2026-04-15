// in order to fit a new bin/dataset/whatever, you need to define a fitting parameter array and a recursive fraction array
// these look like this:
/*
	Double_t par[8][3] - a fitting parameter array
		like-sign background exponent
		unlike-sign background exponent
		1. Gauss mean
		2. Gauss mean
		3. Gauss mean
		1. Gauss width
		2. Gauss width
		3. Gauss width

	Double_t frac[4][3] - a recursive fraction array
		LS background fraction
		US background fraction
		1. Gauss fraction
		2. Gauss fraction
*/
// where in each row you specify: initial fit value; lower constraint; upper constraint

/////////////////////////////////////////////////////////////
// pT and TofMult integrated
vector<vector<Double_t>> parBHT2_int =
{
	{-0.232, -0.5, -0.0},				// LS bg exponent
	{-0.21, -0.5, -0.0},				// US bg exponent
	{9.31, 9.2, 9.4},					// sig1 mean
	{9.9, 9.8, 10.1},					// sig2 mean
	{10.27, 10.2, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_int =
{
	
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
	/*
	2 sigma
	{0.7126, 0.7074, 0.7178},	// sig1_frac
	{0.6828, 0.6782, 0.6874}	// sig2_frac
	*/
};

map_par["BHT2_int"]  = parBHT2_int;
map_frac["BHT2_int"]  = fracBHT2_int;

map_par["BHT2_m1ymp5"]  = parBHT2_int;
map_frac["BHT2_m1ymp5"]  = fracBHT2_int;

map_par["BHT2_mp5yp5"]  = parBHT2_int;
map_frac["BHT2_mp5yp5"]  = fracBHT2_int;

map_par["BHT2_p5y1"]  = parBHT2_int;
map_frac["BHT2_p5y1"]  = fracBHT2_int;

/////////////////////////////////////////////////////////////
// 0 <= TofMult < 3
vector<vector<Double_t>> parBHT2_0mult3 =
{
	{-0.2, -0.5, -0.0},					// LS bg exponent
	{-0.2, -0.5, -0.0},					// US bg exponent
	{9.3, 9.3, 9.35},					// sig1 mean
	{9.9, 9.9, 10.1},					// sig2 mean
	{10.3, 10.2, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_0mult3 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_0mult3"]  = parBHT2_0mult3;
map_frac["BHT2_0mult3"]  = fracBHT2_0mult3;

/////////////////////////////////////////////////////////////
// 3 <= TofMult < 6
vector<vector<Double_t>> parBHT2_3mult6 =
{
	{-0.25, -0.5, -0.0},				// LS bg exponent
	{-0.2, -0.5, -0.0},					// US bg exponent
	{9.33, 9.0, 9.5},					// sig1 mean
	{9.95, 9.9, 10.0},					// sig2 mean
	{10.35, 10.25, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_3mult6 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_3mult6"]  = parBHT2_3mult6;
map_frac["BHT2_3mult6"]  = fracBHT2_3mult6;

/////////////////////////////////////////////////////////////
// 6 <= TofMult < 9
vector<vector<Double_t>> parBHT2_6mult9 =
{
	{-0.2, -0.5, -0.0},					// LS bg exponent
	{-0.2, -0.5, -0.0},					// US bg exponent
	{9.3, 9.0, 9.5},					// sig1 mean
	{9.9, 9.9, 10.0},					// sig2 mean
	{10.3, 10.25, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_6mult9 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_6mult9"]  = parBHT2_6mult9;
map_frac["BHT2_6mult9"]  = fracBHT2_6mult9;

/////////////////////////////////////////////////////////////
// 9 <= TofMult < 12
vector<vector<Double_t>> parBHT2_9mult12 =
{
	{-0.2, -0.5, -0.0},					// LS bg exponent
	{-0.2, -0.5, -0.0},					// US bg exponent
	{9.442, 9.0, 9.5},					// sig1 mean
	{10.006, 9.9, 10.1},				// sig2 mean
	{10.334, 10.2, 10.5},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_9mult12 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_9mult12"]  = parBHT2_9mult12;
map_frac["BHT2_9mult12"]  = fracBHT2_9mult12;

/////////////////////////////////////////////////////////////
// 12 <= TofMult < 20
vector<vector<Double_t>> parBHT2_12mult20 =
{
	{-0.2, -0.5, -0.0},					// LS bg exponent
	{-0.2, -0.5, -0.0},					// US bg exponent
	{9.442, 9.0, 9.5},					// sig1 mean
	{10.006, 9.9, 10.1},				// sig2 mean
	{10.334, 10.2, 10.5},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_12mult20 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_12mult20"]  = parBHT2_12mult20;
map_frac["BHT2_12mult20"]  = fracBHT2_12mult20;

/////////////////////////////////////////////////////////////
// 9 <= TofMult < 20
vector<vector<Double_t>> parBHT2_9mult20 =
{
	{-0.2, -0.5, -0.0},					// LS bg exponent
	{-0.2, -0.5, -0.0},					// US bg exponent
	{9.3, 9.0, 9.5},					// sig1 mean
	{9.9, 9.9, 9.95},					// sig2 mean
	{10.3, 10.25, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_9mult20 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_9mult20"]  = parBHT2_9mult20;
map_frac["BHT2_9mult20"]  = fracBHT2_9mult20;

/////////////////////////////////////////////////////////////
// 0 <= pT < 2
vector<vector<Double_t>> parBHT2_0pT2 =
{
	{-0.232, -0.5, 0.5},				// LS bg exponent
	{-0.21, -0.5, 0.5},				// US bg exponent
	{9.31, 9.2, 9.4},					// sig1 mean
	{9.9, 9.8, 10.1},					// sig2 mean
	{10.27, 10.2, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_0pT2 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_0pT2"]  = parBHT2_0pT2;
map_frac["BHT2_0pT2"]  = fracBHT2_0pT2;

/////////////////////////////////////////////////////////////
// 2 <= pT < 4
vector<vector<Double_t>> parBHT2_2pT4 =
{
	{-0.232, -0.5, -0.0},				// LS bg exponent
	{-0.21, -0.5, -0.0},				// US bg exponent
	{9.31, 9.2, 9.4},					// sig1 mean
	{9.9, 9.8, 10.1},					// sig2 mean
	{10.27, 10.2, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_2pT4 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_2pT4"]  = parBHT2_2pT4;
map_frac["BHT2_2pT4"]  = fracBHT2_2pT4;

/////////////////////////////////////////////////////////////
// 4 <= pT < 6
vector<vector<Double_t>> parBHT2_4pT6 =
{
	{-0.232, -0.5, -0.0},				// LS bg exponent
	{-0.21, -0.5, -0.0},				// US bg exponent
	{9.31, 9.2, 9.4},					// sig1 mean
	{9.9, 9.8, 10.1},					// sig2 mean
	{10.27, 10.2, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_4pT6 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_4pT6"]  = parBHT2_4pT6;
map_frac["BHT2_4pT6"]  = fracBHT2_4pT6;

/////////////////////////////////////////////////////////////
// 6 <= pT < 8
vector<vector<Double_t>> parBHT2_6pT8 =
{
	{-0.232, -0.5, -0.0},				// LS bg exponent
	{-0.21, -0.5, -0.0},				// US bg exponent
	{9.31, 9.2, 9.4},					// sig1 mean
	{9.9, 9.8, 10.1},					// sig2 mean
	{10.27, 10.2, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_6pT8 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_6pT8"]  = parBHT2_6pT8;
map_frac["BHT2_6pT8"]  = fracBHT2_6pT8;

/////////////////////////////////////////////////////////////
// 8 <= pT < 10
vector<vector<Double_t>> parBHT2_8pT10 =
{
	{-0.232, -0.5, -0.0},				// LS bg exponent
	{-0.21, -0.5, -0.0},				// US bg exponent
	{9.31, 9.2, 9.4},					// sig1 mean
	{9.9, 9.8, 10.1},					// sig2 mean
	{10.27, 10.2, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_8pT10 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};

map_par["BHT2_8pT10"]  = parBHT2_8pT10;
map_frac["BHT2_8pT10"]  = fracBHT2_8pT10;
/////////////////////////////////////////////////////////////
// 10 <= pT < 15
vector<vector<Double_t>> parBHT2_10pT15 =
{
	{-0.232, -0.5, -0.0},				// LS bg exponent
	{-0.21, -0.5, -0.0},				// US bg exponent
	{9.31, 9.2, 9.4},					// sig1 mean
	{9.9, 9.8, 10.1},					// sig2 mean
	{10.27, 10.2, 10.35},				// sig3 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n
	{1.678214, 1.615692, 1.740736},		// sig2 n
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> fracBHT2_10pT15 =
{
	{0.74, 0.4, 1.0},		// LS bg frac
    {0.63, 0.0, 1.0},		// US bg frac
    {0.7126, 0.71, 0.715},	// sig1_frac
    {0.6828, 0.68, 0.685}	// sig2_frac
};