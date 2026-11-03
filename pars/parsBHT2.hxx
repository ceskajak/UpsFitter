/////////////////////////////////////////////////////////////
// pT and TofMult integrated
vector<vector<Double_t>> par_sig_BHT2_int =
{
	{9.31, 9.2, 9.6},					// sig1 mean (upper bound widened from 9.4: semi-free, PDG 9.46)
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

vector<vector<Double_t>> par_bg_BHT2_cheby =
{
	{-0.6, -1.0, 1.0},					// comb/LS bg Chebyshev c1
	{0.1, -1.0, 1.0},					// comb/LS bg Chebyshev c2
	{-0.6, -1.0, 1.0},					// US bg Chebyshev c1
	{0.1, -1.0, 1.0}					// US bg Chebyshev c2
};

vector<vector<Double_t>> frac_BHT2_int =
{
	{0.74, 0.0, 1.0},					// comb/LS bg frac (fixed to LS/US at fit time, range only has to contain it)
    {0.63, 0.0, 1.0},					// res bg frac
    {0.7126, 0.71, 0.715},				// sig1 frac
    {0.6828, 0.68, 0.685}				// sig2 frac
};

map_par_sig["BHT2_int"]  = par_sig_BHT2_int;
map_par_bg["BHT2_exp"] = par_bg_BHT2_exp;
map_par_bg["BHT2_cheby"] = par_bg_BHT2_cheby;
map_frac["BHT2_int"]  = frac_BHT2_int;

/////////////////////////////////////////////////////////////
// Multiplicity bins
// Lineshape and fraction start values and bounds are those of BHT2_int (embedding), only the background start values differ per bin
// Start values {comb exponent, res exponent, res frac} from the nominal fits (subcorr), bounds as in BHT2_int
map<TString, vector<Double_t>> start_BHT2 =
{
	{"0mult1", {-0.222, -0.175, 0.472}},
	{"1mult2", {-0.249, -0.244, 0.635}},
	{"2mult3", {-0.240, -0.279, 0.719}},
	{"3mult4", {-0.248, -0.205, 0.737}},
	{"4mult8", {-0.259, -0.166, 0.856}},
	{"3mult8", {-0.253, -0.195, 0.783}}
};

// Chebyshev start values {comb c1, comb c2, res c1, res c2}, bounds as in par_bg_BHT2_cheby
map<TString, vector<Double_t>> start_BHT2_cheby =
{
	{"0mult1", {-0.597, -0.018, -0.461, -0.079}},
	{"1mult2", {-0.633, 0.087, -0.591, 0.1}},
	{"2mult3", {-0.603, 0.132, -0.722, 0.030}},
	{"3mult4", {-0.627, 0.108, -0.535, 0.095}},
	{"4mult8", {-0.660, 0.102, -0.470, -0.288}},
	{"3mult8", {-0.642, 0.106, -0.522, -0.036}}
};

for (auto & start : start_BHT2_cheby)
{
	vector<vector<Double_t>> par_bg = par_bg_BHT2_cheby;
	for (int i = 0; i < 4; i++) par_bg[i][0] = start.second[i];
	map_par_bg["BHT2_cheby_" + start.first] = par_bg;
}

for (auto & start : start_BHT2)
{
	vector<vector<Double_t>> par_bg = par_bg_BHT2_exp;
	vector<vector<Double_t>> frac = frac_BHT2_int;
	par_bg[0][0] = start.second[0];
	par_bg[1][0] = start.second[1];
	frac[1][0] = start.second[2];

	map_par_sig["BHT2_" + start.first] = par_sig_BHT2_int;
	map_par_bg["BHT2_exp_" + start.first] = par_bg;
	map_frac["BHT2_" + start.first] = frac;
}

