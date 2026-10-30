/////////////////////////////////////////////////////////////
// pT and TofMult integrated
vector<vector<Double_t>> par_sig_BHT2_int =
{
	{9.31, 9.2, 9.4},					// sig1 mean
	{0.302521, 0.300231, 0.304811},		// sig1 width
	{1.639359, 1.610738, 2.183158},		// sig1 alpha
	{1.714986, 1.644381, 1.785591},		// sig1 n

	{9.9, 9.8, 10.1},					// sig2 mean
	{0.323332, 0.321315, 0.325349},		// sig2 width
	{1.649815, 1.622159, 1.677471},		// sig2 alpha
	{1.678214, 1.615692, 1.740736},		// sig2 n

	{10.27, 10.2, 10.35},				// sig3 mean
	{0.337254, 0.335369, 0.339139},		// sig3 width
	{1.653388, 1.630794, 1.675982},		// sig3 alpha
	{1.657812, 1.606731, 1.708893}		// sig3 n
};

vector<vector<Double_t>> par_bg_BHT2_exp =
{
	{-0.232, -0.5, -0.0},				// comb/LS bg exponent
	{-0.21, -0.5, -0.0}					// US bg exponent
};

vector<vector<Double_t>> frac_BHT2_int =
{
	{0.74, 0.4, 1.0},					// comb/LS bg frac
    {0.63, 0.0, 1.0},					// res bg frac
    {0.7126, 0.71, 0.715},				// sig1 frac
    {0.6828, 0.68, 0.685}				// sig2 frac
};

map_par_sig["BHT2_int"]  = par_sig_BHT2_int;
map_par_bg["BHT2_exp"] = par_bg_BHT2_exp;
map_frac["BHT2_int"]  = frac_BHT2_int;

/////////////////////////////////////////////////////////////