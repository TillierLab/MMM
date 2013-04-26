#include "DistanceMatrix.hpp"

double DistanceMatrix::compute_pmb_distance_quick(vector <unsigned char> const& sequence1, vector <unsigned char> const& sequence2)
{	double distance;
	bool overlap = false;

	unsigned int total=0, different=0, length=sequence1.size();
	for (unsigned int pos = 0; pos < length; pos++)
	{	unsigned char const& aa1 = sequence1.at(pos);
		unsigned char const& aa2 = sequence2.at(pos);
		
		//skip non-aminoacid positions
		if (aa1 == SKIP_AA || aa2 == SKIP_AA)
			continue;
		overlap=true;
		total++;
		if(aa1 != aa2)
			different++;
	}
	if (!overlap)
	{	//printf("\nWARNING: NO OVERLAP BETWEEN SEQUENCES %ld AND %ld; -1.0 WAS WRITTEN\n", 0, 0);
		distance = DISTANCE_NO_OVERLAP;
	} else
	{	double distanceObserved = different/(double)total;
		if(distanceObserved < 0.8813)
		{	distance = (1.228 * distanceObserved) / (pow((1.744+distanceObserved),0.365) * pow((1.004-distanceObserved), 0.635)); //Formula from "Veerassamy et al (2003) A transition Probability model for amino acid substitutions from blocks"
		} else
		{	//printf("\nWARNING: THE DISTANCE EXCEEDS THE THRESHOLD FOR THE SEQUENCES %ld AND %ld; 0.0 WAS WRITTEN\n", 0, 0);
			distance = DISTANCE_INFINITY;
		}
	}
	return distance;
}
//


double DistanceMatrix::compute_pmb_distance(vector <unsigned char> const& sequence1, vector <unsigned char> const& sequence2)
{	static const int N_OF_AAS = 20;
	
		/* PMB matrix decomposition courtesy of Elisabeth Tillier */
	static const double eig[N_OF_AAS] = 
	{0.0000001586972220,-1.8416770496147100, -1.6025046986139100,-1.5801012515121300,
	-1.4987794099715900,-1.3520794233801900,-1.3003469390479700,-1.2439503327631300,
	-1.1962574080244200,-1.1383730501367500,-1.1153278910708000,-0.4934843510654760,
	-0.5419014550215590,-0.9657997830826700,-0.6276075673757390,-0.6675927795018510,
	-0.6932641383465870,-0.8897872681859630,-0.8382698977371710,-0.8074694642446040};
	//
	
	static const double prob[N_OF_AAS][N_OF_AAS] =
	{{0.0771762457248147,0.0531913844998640,0.0393445076407294,0.0466756566755510,
	0.0286348361997465,0.0312327748383639,0.0505410248721427,0.0767106611472993,
	0.0258916271688597,0.0673140562194124,0.0965705469252199,0.0515979465932174,
	0.0250628079438675,0.0503492018628350,0.0399908189418273,0.0641898881894471,
	0.0517539616710987,0.0143507440546115,0.0357994592438322,0.0736218495862984},
	{0.0368263046116572,-0.0006728917107827,0.0008590805287740,-0.0002764255356960,
	0.0020152937187455,0.0055743720652960,0.0003213317669367,0.0000449190281568,
	-0.0004226254397134,0.1805040629634510,-0.0272246813586204,0.0005904606533477,
	-0.0183743200073889,-0.0009194625608688,0.0008173657533167,-0.0262629806302238,
	0.0265738757209787,0.0002176606241904,0.0021315644838566,-0.1823229927207580},
	{-0.0194800075560895,0.0012068088610652,-0.0008803318319596,-0.0016044273960017,
	-0.0002938633803197,-0.0535796754602196,0.0155163896648621,-0.0015006360762140,
	0.0021601372013703,0.0268513218744797,-0.1085292493742730,0.0149753083138452,
	0.1346457366717310,-0.0009371698759829,0.0013501708044116,0.0346352293103622,
	-0.0276963770242276,0.0003643142783940,0.0002074817333067,-0.0174108903914110},
	{0.0557839400850153,0.0023271577185437,0.0183481103396687,0.0023339480096311,
	0.0002013267015151,-0.0227406863569852,0.0098644845475047,0.0064721276774396,
	0.0001389408104210,-0.0473713878768274,-0.0086984445005797,0.0026913674934634,
	0.0283724052562196,0.0001063665179457,0.0027442574779383,-0.1875312134708470,
	0.1279864877057640,0.0005103347834563,0.0003155113168637,0.0081451082759554},
	{0.0037510125027265,0.0107095920636885,0.0147305410328404,-0.0112351252180332,
	-0.0001500408626446,-0.1523450933729730,0.0611532413339872,-0.0005496748939503,
	0.0048714378736644,-0.0003826320053999,0.0552010244407311,0.0482555671001955,
	-0.0461664995115847,-0.0021165008617978,-0.0004574454232187,0.0233755883688949,
	-0.0035484915422384,0.0009090698422851,0.0013840637687758,-0.0073895139302231},
	{-0.0111512564930024,0.1025460064723080,0.0396772456883791,-0.0298408501361294,
	-0.0001656742634733,-0.0079876311843289,0.0712644184507945,-0.0010780604625230,
	-0.0035880882043592,0.0021070399334252,0.0016716329894279,-0.1810123023850110,
	0.0015141703608724,-0.0032700852781804,0.0035503782441679,0.0118634302028026,
	0.0044561606458028,-0.0001576678495964,0.0023470722225751,-0.0027457045397157},
	{0.1474525743949170,-0.0054432538500293,0.0853848892349828,-0.0137787746207348,
	-0.0008274830358513,0.0042248844582553,0.0019556229305563,-0.0164191435175148,
	-0.0024501858854849,0.0120908948084233,-0.0381456105972653,0.0101271614855119,
	-0.0061945941321859,0.0178841099895867,-0.0014577779202600,-0.0752120602555032,
	-0.1426985695849920,0.0002862275078983,-0.0081191734261838,0.0313401149422531},
	{0.0542034611735289,-0.0078763926211829,0.0060433542506096,0.0033396210615510,
	0.0013965072374079,0.0067798903832256,-0.0135291136622509,-0.0089982442731848,
	-0.0056744537593887,-0.0766524225176246,0.1881210263933930,-0.0065875518675173,
	0.0416627569300375,-0.0953804133524747,-0.0012559228448735,0.0101622644292547,
	-0.0304742453119050,0.0011702318499737,0.0454733434783982,-0.1119239362388150},
	{0.1069409037912470,0.0805064400880297,-0.1127352030714600,0.1001181253523260,
	-0.0021480427488769,-0.0332884841459003,-0.0679837575848452,-0.0043812841356657,
	0.0153418716846395,-0.0079441315103188,-0.0121766182046363,-0.0381127991037620,
	-0.0036338726532673,0.0195324059593791,-0.0020165963699984,-0.0061222685010268,
	-0.0253761448771437,-0.0005246410999057,-0.0112205170502433,0.0052248485517237},
	{-0.0325247648326262,0.0238753651653669,0.0203684886605797,0.0295666232678825,
	-0.0003946714764213,-0.0157242718469554,-0.0511737848084862,0.0084725632040180,
	-0.0167068828528921,0.0686962159427527,-0.0659702890616198,-0.0014289912494271,
	-0.0167000964093416,-0.1276689083678200,0.0036575057830967,-0.0205958145531018,
	0.0000368919612829,0.0014413626622426,0.1064360941926030,0.0863372661517408},
	{-0.0463777468104402,0.0394712148670596,0.1118686750747160,0.0440711686389031,
	-0.0026076286506751,-0.0268454015202516,-0.1464943067133240,-0.0137514051835380,
	-0.0094395514284145,-0.0144124844774228,0.0249103379323744,-0.0071832157138676,
	0.0035592787728526,0.0415627419826693,0.0027040097365669,0.0337523666612066,
	0.0316121324137152,-0.0011350177559026,-0.0349998884574440,-0.0302651879823361},
	{0.0142360925194728,0.0413145623127025,0.0324976427846929,0.0580930922002398,
	-0.0586974207121084,0.0202001168873069,0.0492204086749069,0.1126593173463060,
	0.0116620013776662,-0.0780333711712066,-0.1109786767320410,0.0407775100936731,
	-0.0205013161312652,-0.0653458585025237,0.0347351829703865,0.0304448983224773,
	0.0068813748197884,-0.0189002309261882,-0.0334507528405279,-0.0668143558699485},
	{-0.0131548829657936,0.0044244322828034,-0.0050639951827271,-0.0038668197633889,
	-0.1536822386530220,0.0026336969165336,0.0021585651200470,-0.0459233839062969,
	0.0046854727140565,0.0393815434593599,0.0619554007991097,0.0027456299925622,
	0.0117574347936383,0.0373018612990383,0.0024818527553328,-0.0133956606027299,
	-0.0020457128424105,0.0154178819990401,0.0246524142683911,0.0275363065682921},
	{-0.1542307272455030,0.0364861558267547,-0.0090880407008181,0.0531673937889863,
	0.0157585615170580,0.0029986538457297,0.0180194047699875,0.0652152443589317,
	0.0266842840376180,0.0388457366405908,0.0856237634510719,0.0126955778952183,
	0.0099593861698250,-0.0013941794862563,0.0294065511237513,-0.1151906949298290,
	-0.0852991447389655,0.0028699120202636,-0.0332087026659522,0.0006811857297899},
	{0.0281300736924501,-0.0584072081898638,-0.0178386569847853,-0.0536470338171487,
	-0.0186881656029960,-0.0240008730656106,-0.0541064820498883,0.2217137098936020,
	-0.0260500001542033,0.0234505236798375,0.0311127151218573,-0.0494139126682672,
	0.0057093465049849,0.0124937286655911,-0.0298322975915689,0.0006520211333102,
	-0.0061018680727128,-0.0007081999479528,-0.0060523759094034,0.0215845995364623},
	{0.0295321046399105,-0.0088296411830544,-0.0065057049917325,-0.0053478115612781,
	-0.0100646496794634,-0.0015473619084872,0.0008539960632865,-0.0376381933046211,
	-0.0328135588935604,0.0672161874239480,0.0667626853916552,-0.0026511651464901,
	0.0140451514222062,-0.0544836996133137,0.0427485157912094,0.0097455780205802,
	0.0177309072915667,-0.0828759701187452,-0.0729504795471370,0.0670731961252313},
	{0.0082646581043963,-0.0319918630534466,-0.0188454445200422,-0.0374976353856606,
	0.0037131290686848,-0.0132507796987883,-0.0306958830735725,-0.0044119395527308,
	-0.0140786756619672,-0.0180512599925078,-0.0208243802903953,-0.0232202769398931,
	-0.0063135878270273,0.0110442171178168,0.1824538048228460,-0.0006644614422758,
	-0.0069909097436659,0.0255407650654681,0.0099119399501151,-0.0140911517070698},
	{0.0261344441524861,-0.0714454044548650,0.0159436926233439,0.0028462736216688,
	-0.0044572637889080,-0.0089474834434532,-0.0177570282144517,-0.0153693244094452,
	0.1160919467206400,0.0304911481385036,0.0047047513411774,-0.0456535116423972,
	0.0004491494948617,-0.0767108879444462,-0.0012688533741441,0.0192445965934123,
	0.0202321954782039,0.0281039933233607,-0.0590403018490048,0.0364080426546883},
	{0.0115826306265004,0.1340228176509380,-0.0236200652949049,-0.1284484655137340,
	-0.0004742338006503,0.0127617346949511,-0.0428560878860394,0.0060030732454125,
	0.0089182609926781,0.0085353834972860,0.0048464809638033,0.0709740071429510,
	0.0029940462557054,-0.0483434904493132,-0.0071713680727884,-0.0036840391887209,
	0.0031454003250096,0.0246243550241551,-0.0449551277644180,0.0111449232769393},
	{0.0140356721886765,-0.0196518236826680,0.0030517022326582,0.0582672093364850,
	-0.0000973895685457,0.0021704767224292,0.0341806268602705,-0.0152035987563018,
	-0.0903198657739177,0.0259623214586925,0.0155832497882743,-0.0040543568451651,
	0.0036477631918247,-0.0532892744763217,-0.0142569373662724,0.0104500681408622,
	0.0103483945857315,0.0679534422398752,-0.0768068882938636,0.0280289727046158}}
	;
	//

	
	double distance = 0.1;
	double delta = distance / 2.0;
	
	unsigned int length = sequence1.size();//possible that it takes a long time? No, doesn't make a difference
	for( unsigned char iterations = 0;  iterations < 20; iterations++)
	{	double slope = 0.0;
		double curv = 0.0;
		bool neginfinity = false;
		bool overlap = false;

		for (unsigned int pos = 0; pos < length; pos++)
		{	unsigned char const& aa1 = sequence1.at(pos);
			unsigned char const& aa2 = sequence2.at(pos);

			//skip non-aminoacid positions
			if (aa1 == SKIP_AA || aa2 == SKIP_AA)
				continue;
			overlap = true;
		
		
			double p = 0.0;
			double dp = 0.0;
			double d2p = 0.0;
		 	for (unsigned char aa_index = 0; aa_index < N_OF_AAS; aa_index++)
			{	double const& eigen = eig[aa_index];
				
				if (gamma_)
				{	double corr = 1.0-eigen*distance/alpha_;
					double corr_eigen = eigen/corr;
					double tmp = prob[aa_index][aa1] * prob[aa_index][aa2] * exp(-alpha_*log(corr));
					
					p += tmp;
					
					tmp*=corr_eigen;
					dp += tmp; 
					
					tmp*=inv_alpha_*corr_eigen;
					d2p += tmp;
					
				} else
				{
					double tmp = prob[aa_index][aa1] * prob[aa_index][aa2] * exp(eigen*distance);
						
					p += tmp;
					
					tmp*=eigen;
					dp += tmp;
					
					tmp*=eigen;
					d2p += tmp;
				}
			}
		
			//printf("pos=%ld\tp=%f\n", pos, p);		
			if (p <= 0.0)
			{	neginfinity = true; 
				break; //should abort the iteration immediately? 
			}
			else
			{	
				
				//double const& dp_p = dp / p;
				double dp_p = dp / p;
				slope += dp_p;
				curv += (d2p / p - dp_p*dp_p);
			}
		}
	 
		if (!overlap)
		{	//printf("\nWARNING: NO OVERLAP BETWEEN SEQUENCES %ld AND %ld; -1.0 WAS WRITTEN\n", i+1, j+1);
			distance = DISTANCE_NO_OVERLAP;
			break;
		} else if (!neginfinity)
		{	if (curv < 0.0)
			{	distance -= slope / curv;
				if (distance > 10000.0)
				{	//printf("\nWARNING: INFINITE DISTANCE BETWEEN SPECIES %ld AND %ld; -1.0 WAS WRITTEN\n", i+1, j+1);
					distance = DISTANCE_INFINITY;
					break;
				}
			}
			else
			{	if ((slope > 0.0 && delta < 0.0) || (slope < 0.0 && delta > 0.0))
					delta /= -2;
				distance += delta;
			}
		} else
		{	delta /= -2;
			distance += delta;
		}
		if (distance < DISTANCE_TOO_SMALL)
			distance =  DISTANCE_TOO_SMALL;
	}
	return distance;
}
//

void DistanceMatrix::create_from_alignment(Alignment alignment, bool quick)
{	static const unsigned char aa2index[26] = {ala, SKIP_AA, cys, asp, glu, phe, gly, his, ileu, SKIP_AA, lys, leu, met, asn, SKIP_AA, pro, gln, arg, ser, thr, SKIP_AA, val, trp, SKIP_AA, tyr, SKIP_AA};
	
	init();
		
	nOfSequences_ = alignment.get_nOfSequences();
	int seqLength=alignment.get_nOfColumns();
	headers_.reserve(nOfSequences_);
	distances2D_.resize(nOfSequences_);
		
	//calcualte pairwise distances for all sequences
	vector<vector<unsigned char > > sequences;
	sequences.resize(nOfSequences_);
	for (unsigned int i=0; i<nOfSequences_;i++)
	{	headers_.push_back(alignment.get_header(i));
		
		distances2D_.at(i).resize(nOfSequences_);
		distances2D_.at(i).at(i)=0.0;
				
		
		string const& sequenceStr = alignment.get_sequence(i);
		vector<unsigned char>& sequence1 = sequences.at(i);
		sequence1.reserve(seqLength);
		//convert aminoacid codes to integers (for direct correspondence with scoring matrices)
		for (int pos=0; pos<seqLength; pos++)
		{	sequence1.push_back((isupper(sequenceStr.at(pos)))? aa2index[sequenceStr.at(pos)-'A'] : SKIP_AA);
		}
		
		
		//compute pairwise distances for this sequence
		for (unsigned int j=0; j<i;j++)
		{	vector <unsigned char> const&  sequence2 = sequences.at(j);
			double distance = (quick)? compute_pmb_distance_quick(sequence1, sequence2) : compute_pmb_distance(sequence1, sequence2);
			distances2D_.at(i).at(j) = distance;
			distances2D_.at(j).at(i) = distance;
		}
	}
	valid_ = true;
}
//

void DistanceMatrix::read_from_filehandle(ifstream& inFile)
{	init();
	if (!inFile.is_open())
	{	throw_error("invalid filehandle");
	}
	
	inFile >> nOfSequences_;
	if (inFile.peek() != '\n')
		throw_error("invalid distance matrix format 1");
	
	headers_.reserve(nOfSequences_);
	distances2D_.resize(nOfSequences_);
	
	for (unsigned int i = 0; i < nOfSequences_; i++)
	{	string header;
		inFile >> header;
		headers_.push_back(header);
		//taxons_.push_back(header.substr(0, header.find('|'))); // If not found, selects the whole name
		
		vector <double>& distancesRow = distances2D_.at(i);
		distancesRow.reserve(nOfSequences_);
		for (unsigned int j = 0; j < nOfSequences_; j++)
		{	double distance;
			inFile >> distance;
			distancesRow.push_back(max(1e-20,distance)); // Don't allow negative or zero distances
		}
		char next = inFile.get();
		if (next != '\n' || distancesRow.size() != nOfSequences_)
		{	//throw runtime_error(get_name() + "invalid distance matrix format 2");
			stringstream err;
			err << " Read " << distancesRow.size() << "/" << nOfSequences_ << " distances. Next: [" << next <<"]";
			throw_error("invalid distance matrix format 2" + err.str());
		}
	}
	if (headers_.size() != nOfSequences_ || !inFile.good())
			throw_error("invalid distance matrix format 3");
	inFile.peek();
	if (!inFile.eof())
		throw_error("invalid distance matrix format 4");
	valid_ = true;
}
//

void DistanceMatrix::read_from_file(string const& inFilepath)
{	init();
	ifstream inFile(inFilepath.c_str());
	if (!inFile.is_open())
	{	//cout << "\nError: cannot open file " << inFilepath << endl;
		//exit(1);
		throw_error("cannot open file: " + inFilepath);
	}
	
	read_from_filehandle(inFile);
	inFile.close();
}
//

void DistanceMatrix::print_to_filehandle(ofstream& outFile) const
{	validate();
	if (!outFile.is_open())
	{	//cout << "\nError: invalid filehandle" << endl;
		//exit(1);
		throw_error("invalid filehandle");
	}
	
	//print the number of sequences
	outFile.width(5); 
	outFile << nOfSequences_ << '\n';
		
	outFile.setf(ios::fixed,ios::floatfield);
	outFile.precision(6);
	//output distances
	for (unsigned int i = 0; i < nOfSequences_; i++)
	{	//print the name
		outFile << headers_.at(i).substr(0,10);
		for (unsigned int j = 0; j < nOfSequences_; j++)
		{	double const& dst = distances2D_.at(i).at(j);
			if (dst < 10.0)
				outFile << ' ';
			outFile << ' ' << dst;
			
			
			//	fprintf(outfile, "%10.6f", dst);
			//else if (dst < 1000.0)
			//	fprintf(outfile, " %10.6f", dst);
			//else 
			//	fprintf(outfile, " %11.6f", dst);
			
			if ((j + 2) % 7 == 0 && j+1 < nOfSequences_)
				outFile << '\n';
		}
		outFile << '\n';
	}
}
//

void DistanceMatrix::print_to_file(string const& outFilepath) const
{	validate();
	ofstream outFile(outFilepath.c_str());
	if (!outFile.is_open()) {
		//cout << "\nError: cannot open file " << outFilepath << endl;
		//exit(1);
		throw_error("cannot create file: " + outFilepath);
	}
	print_to_filehandle(outFile);
	outFile.close();
}
//

void DistanceMatrix::get_distances_as_1d(vector <double>& distances1D) const
{	validate();
	distances1D.clear();
	distances1D.reserve(nOfSequences_*nOfSequences_);
	for (unsigned int i=0; i<nOfSequences_; i++)
	{	for (unsigned int j=0; j<nOfSequences_; j++)
		{	double dstRounded;
			
			//reduce precision for compatibility
			{	stringstream tmp;
				tmp.setf(ios::fixed,ios::floatfield);
				tmp.precision(6);
				tmp << distances2D_.at(i).at(j);
				tmp >> dstRounded;
			}
			distances1D.push_back(dstRounded);
			
			
			//distances1D.push_back(distances2D_.at(i).at(j));
		}
	}
	return;
}

//
