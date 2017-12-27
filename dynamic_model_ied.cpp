#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <algorithm>
#include <ctime>
#include <cmath>
//#include <nlopt.h>
#include <nlopt.hpp>


using namespace std;

typedef vector<vector<double>> d_vec;

//cell unit: k cells
struct parameters
{
	double u; //attachment rate
	double v; //detachment rate
	double min_n; //minimum number of cells to be considered non-zero
	double p0;
	double beta;
	double w;
	double K;
	
	double x;
	double r;
	d_vec S;
	
	int tlag;
	double thalf;
	double Amax;
	double n;

	double phi;
	double psai;
	double theta;
	double C;
	double ecstart;

};

class BB //state of a subpopulation of Bb
{
public:
	int extinct_bool; //if 1, then go extinct this round
	double n_free; //free cell number
	double n_bound; //bound cell number
	int ancestor_type; //type of ancestor vlsE, the root is -1 ***type IS index
	vector<int> nr; //number of recombinations away
};

vector<vector<int>> nr_matrix; //nr matrix

class AB //state of a subpopulation of Ab
{
public:
	int type; //specific vlse type
	double concentration; //concentration of the subpopulation
	int tinit; //initiation time
};

class EC //state of immune effector cells
{
public:
	double concentration;
	
};

int nr(int type1, int type2, vector<BB> bb) //the number of recombinations till common ancestor
{

	int n;

	vector<int> v1;
	vector<int> v2;

	for(int i=type1; bb[i].ancestor_type >=0; i=bb[i].ancestor_type)	v1.insert(v1.begin(),i);
	for(int i=type2; bb[i].ancestor_type >=0; i=bb[i].ancestor_type)	v2.insert(v2.begin(),i);

	int i;
	for(i=0;i<v1.size()&&i<v2.size();i++)
	{
		if(v1[i]!=v2[i]) {break;}
	}

	n=v1.size()+v2.size()-2*i;

return n;
}

double rxn_AB(int type, vector<BB> bb , vector<AB> ab, parameters p) //calculate effective ab population targeting bb type
//***about S vector: it stores 1000 similarity values for each possible nr; in S[0-999][nr-1], but nr-1 in 0-19
{
	double rab=0;
	for(int i=0;i<ab.size();i++)
	{
		int recombination = nr_matrix[type][i]; 
		recombination = ( (recombination>19) ? 19 : recombination );
		int sample = rand()%1000;
		rab=rab+ pow(p.S[sample][recombination],p.x) * ab[i].concentration;
		//cout<<p.S[sample][recombination]<<" "<<(double)p.x<<" "<<ab[i].concentration<<endl;
	}
	
return rab;
}

void grow_AB(vector<AB> &ab, vector<BB> bb, int t, parameters p) //growth of ab, bb is the snapshot
{
	//add new ab types
	for(int bi=ab.size();bi<bb.size();bi++)
	{
		AB newab;
		newab.type=bi;
		newab.concentration=0.000000001; //may need to adjust..
		newab.tinit=t;
		ab.push_back(newab);
	}

	//calculate unnormalized ab concentration
	for(int ai=0;ai<ab.size();ai++)
	{
		if(t-ab[ai].tinit <= p.tlag) ab[ai].concentration=0.000000001;
		else ab[ai].concentration= p.Amax * pow((t-ab[ai].tinit-p.tlag),p.n) / (pow(p.thalf,p.n) + pow((t-ab[ai].tinit-p.tlag),p.n));
	}

}

double Aob(int t, parameters p) //return the observed value of Ab density at time t	
{
	double aob;
	if(t <= p.tlag) aob=0.0;
	else aob=p.Amax * pow((t-p.tlag),p.n) / (pow(p.thalf,p.n) + pow((t-p.tlag),p.n));

return aob;
}

void normalize_AB(vector<AB> &ab, int t, parameters p) //normalize based on Ab ceiling
{
	double abtotal=0.0;

	for(int ai=0;ai<ab.size();ai++)
	{	
		abtotal=abtotal+ab[ai].concentration;
	}

	if(abtotal!=0)
	{
		for(int ai=0;ai<ab.size();ai++) ab[ai].concentration=ab[ai].concentration*Aob(t,p)/abtotal;
	}

}

double get_total_ab_size(vector<AB> ab)
{
	double abtotal=0.0;

	for(int ai=0;ai<ab.size();ai++)
	{	
		abtotal=abtotal+ab[ai].concentration;
	}
return abtotal;
}

double get_total_bb_size(vector<BB> bb)
{
	double bb_total=0.0;
	for(int i=0;i<bb.size();++i)
	{
		if(bb[i].extinct_bool==0) 
		{
			bb_total+=bb[i].n_free+bb[i].n_bound;
		}
	}
	
return bb_total;
}


void bind_grow_recombine_BB(vector<BB> &bb, EC ec , vector<AB> ab, parameters p, int t) //binding+growing+recombining bb
{
	vector<BB> present_bb=bb; //snapshot

	//================calculate logistic growth rate====================
	double total_n_free=0;
	for(int i=0;i<bb.size();i++) total_n_free=total_n_free+present_bb[i].n_free;
	double growth_rate=p.w*p.beta*(1.0-total_n_free/p.K);
	
	//==================================================================
	
	//new bb types from recombination
	vector<BB> new_bb;

	//old bb types:


	double total_recombine=0.0; //If too slow, only allow to add one new type every hour.

	for(int i=0;i<bb.size();i++)
	{
		if(bb[i].extinct_bool==1) continue;


		
		
		double recombine_b = 0.0;
		double recombine_f = 0.0;
		double recombine=0.0;
		if(t%24==0) // recombine every six hours
		{
			recombine_b = p.r * present_bb[i].n_bound; 
			recombine_f = p.r * present_bb[i].n_free;
			recombine = recombine_b+recombine_f;
			if(recombine <= p.min_n) recombine=0.0;

			
		}


		//add new bb from recombination

		if(recombine!=0.0)
		{
			BB newbb;
			newbb.n_free=recombine;
			newbb.n_bound=0;
			newbb.extinct_bool=0;
			newbb.ancestor_type=i;

			//update nr
			vector<int> newnr;
			for(int ni=0;ni<nr_matrix.size();ni++)
			{
				//parental distance+1:
				int cd=nr_matrix[ni][i]+1;
				nr_matrix[ni].push_back(cd);
				newnr.push_back(cd);
			}
			newnr.push_back(0);
			nr_matrix.push_back(newnr);
			
			new_bb.push_back(newbb);
		}

		//growth+bind+recombine -> death -> detach; because death and detach can affect the full bound population

		//growth
		double growth = growth_rate * present_bb[i].n_free;
		
		//bind
		double bind = p.u * present_bb[i].n_free * rxn_AB(i,bb,ab,p); //if(bind <= p.min_n) bind=0;
	
		//free cells:
		if(bb[i].extinct_bool!=1)
		bb[i].n_free=bb[i].n_free - bind - recombine_f + growth;

		//bound cells: 
		if(bb[i].extinct_bool!=1)
		bb[i].n_bound=bb[i].n_bound - recombine_b + bind;


		//death
		double death = p.p0 * bb[i].n_bound * ec.concentration ;
		bb[i].n_bound=bb[i].n_bound - death;

		//detach:
		double detach = p.v * present_bb[i].n_bound ; if(detach <= p.min_n) detach=0;
		bb[i].n_free=bb[i].n_free + detach;
		bb[i].n_bound=bb[i].n_bound - detach;
		

		//mark with population extinction:
		if(bb[i].n_free <= p.min_n) bb[i].n_free=0;
		if(bb[i].n_bound <= p.min_n) bb[i].n_bound=0;
		if(bb[i].n_free==0 && bb[i].n_bound==0) bb[i].extinct_bool=1;

	}

	//combine new and old populations
	bb.insert(bb.end(), new_bb.begin(), new_bb.end());


	
	
}
	


void grow_EC(EC &ec, vector<BB> bb, parameters p) //growth of ec based on states of bb
{
	EC present_ec;
	present_ec.concentration=ec.concentration; //snapshot
	
	//total Bb size
	double n_total=0.0;
	for(int bi=0;bi<bb.size();bi++)
	{
		n_total=bb[bi].n_free+bb[bi].n_bound+n_total;
	}

	ec.concentration=ec.concentration + p.phi + p.psai * present_ec.concentration * n_total / (p.C + n_total) - p.theta * present_ec.concentration;
}



double get_nonparental_bb_ratio(vector<BB> bb)
{
	double bb_total=get_total_bb_size(bb);
	double nonparental_bb=0.0;
	double ratio=0.0;

	if(bb_total!=0)
	{
		for(int i=0;i<bb.size();++i)
		{
			if(bb[i].extinct_bool==0&&bb[i].ancestor_type!=-1) 
			{
				nonparental_bb=nonparental_bb+bb[i].n_free+bb[i].n_bound;
			}
		}
		ratio=nonparental_bb/bb_total;
	}
return ratio;
}

void cout_stats(vector<BB> bb, vector<AB> ab, parameters p) //output statistics
{
	int bb_type=0;
	double bb_total=0.0;
	double ab_total=0.0;
	double ab_pressure=0.0;
	
	for(int i=0;i<bb.size();++i)
	{
		if(bb[i].extinct_bool==0) 
		{
			bb_type++;
			bb_total+=bb[i].n_free+bb[i].n_bound;
			ab_pressure=ab_pressure+(bb[i].n_free+bb[i].n_bound)*rxn_AB(i, bb , ab, p);
		}
	}

	for(int i=0;i<ab.size();++i)
	{
		ab_total=ab_total+ab[i].concentration;
	}

	if(bb_total!=0) ab_pressure=ab_pressure/bb_total;
	cout<<bb.size()<<"\t"<<bb_type<<"\t"<<bb_total<<"\t"<<ab_total<<"\t"<<ab_pressure<<endl;
}

typedef struct 
{
	double n_24;
	double n_72;
	double n_192;
	double n_360;
	double n_672;

	double r_96;
	double r_168;
	double r_240;
	double r_336;
	double r_672;

} empirical_data;

double objective(const vector<double> &present_p, vector<double> &grad, void *data) 
{
	empirical_data *d = reinterpret_cast<empirical_data*>(data);

	//========================1. parameter initialization========================
	parameters p;

	p.beta=0.06; //hours -1
	p.K=150; //k cells
	//p.p0=0.02; //cells -1 hours -1
	p.Amax=8.5; //ug ml-1
	p.n=4.0;
	p.u=0.0085; //ml ug-1 hours -1
	p.v=1; //hours
	p.tlag=72; //hours
	p.phi=0.1; // cells hours -1
	p.theta=0.0076; // hours -1
	p.min_n=0.001;
	//p.r=0.0059; //hour-1
	p.r=0.133; //24hours -1
	p.thalf=100.0; //hour
	
	

	p.w=1.0;
	//present_p[0]; //hours -1

	p.psai=present_p[0]; // hours-1
	p.phi=present_p[1]; //cells	
	p.p0=present_p[2]; //cells -1 hours -1

	p.C=present_p[3]; //cells
	p.x=present_p[4];
 


	ifstream S_matrix("S_matrix.txt");
	string line;
	for(;getline(S_matrix,line);)
	{
		vector<double> row_sample;
		string field;
		istringstream linestream(line);
		for(;getline(linestream,field,'\t');) row_sample.push_back(stod(field));
		
		p.S.push_back(row_sample);
	}

	//======2. simulate to capture population size and non-wt vlse type datapoint======



	vector<BB> bb;
	vector<AB> ab;
	EC ec;
	
 	//initialization (hour 0)
	BB newbb;
	newbb.extinct_bool=0;
	newbb.n_free=0.1;
	newbb.n_bound=0;
	newbb.ancestor_type=-1;

	bb.push_back(newbb);
	ec.concentration=0.0;

	vector<int> newnr;
	newnr.push_back(0);
	nr_matrix.push_back(newnr);
	

	//simulation
	double n_24,n_72,n_192,n_360,n_672;
	double r_96,r_168,r_240,r_336,r_672;

	cout<<"=======Simulation========"<<endl;
	cout<<p.psai<<"\t"<<p.phi<<"\t"<<p.p0<<"\t"<<p.C<<"\t"<<p.x<<endl; // hours-1

	for(int hour=1;hour<=720;hour++)
	{
		//make snapshot of all
		vector<BB> present_bb=bb;
		vector<AB> present_ab=ab;
		EC present_ec=ec;

		//grow AB, BB and ec
		grow_AB(ab, present_bb, hour, p);
		normalize_AB(ab, hour, p);
		bind_grow_recombine_BB(bb, present_ec , present_ab, p, hour);
		
		grow_EC(ec, bb, p);

		//record predicted datapoint values

		if(hour==24) n_24=get_total_bb_size(bb);
		if(hour==72) n_72=get_total_bb_size(bb);
		if(hour==192) n_192=get_total_bb_size(bb);
		if(hour==360) n_360=get_total_bb_size(bb);
		//if(hour==672) n_672=get_total_bb_size(bb);

		if(hour==96) r_96=get_nonparental_bb_ratio(bb);
		if(hour==168) r_168=get_nonparental_bb_ratio(bb);
		if(hour==240) r_240=get_nonparental_bb_ratio(bb);
		if(hour==336) r_336=get_nonparental_bb_ratio(bb);
		//if(hour==672) r_672=get_nonparental_bb_ratio(bb);
		if(hour%24==0) cout<<get_total_bb_size(bb)<<"\t"<<get_nonparental_bb_ratio(bb)<<endl;

	}

	//=====================3. compute least square========================
	double square=pow((n_24-d->n_24),2.0)+pow((n_72-d->n_72),2.0)+pow((n_192-d->n_192),2.0)+pow((n_360-d->n_360),2.0)+pow((r_96-d->r_96),2.0)+pow((r_168-d->r_168),2.0)+pow((r_240-d->r_240),2.0)+pow((r_336-d->r_336),2.0);
	cout<<"SD="<<square<<endl;
return square;
}
	


int main()
{

 //=======fitting===========

	//----data-----
	empirical_data data;

	data.n_24=0.5258;
	data.n_72=0.88;
	data.n_192=87;
	data.n_360=6.040;

	data.n_192=87;
	data.n_360=(1.736+6.040)/2;
	data.n_672=10.253;

	data.r_96=0.49;
	data.r_168=0.49;
	data.r_240=0.78;
	data.r_336=0.97;
	data.r_672=0;

	// data source 2:
	// n_192=501.19
	// n_360=10;
	// n_0=10^1.25
	//necessary objects

	//-------------

	nlopt::opt opt(nlopt::LN_COBYLA, 5); /* algorithm and dimensionality */

	opt.set_min_objective(objective, &data); 

	vector<double> present_p(5);
	//present_p[0]=0.9;
	present_p[0]=0.040205; //psai
	present_p[1]=0.516659; //phi
	present_p[2]=0.0159852; //p0

	present_p[3]=1.04632; //C
	present_p[4]=2.55788; //x
	double minf;

	nlopt::result result = opt.optimize(present_p, minf);

	cout<<present_p[0]<<"\t"<<present_p[1]<<"\t"<<present_p[2]<<"\t"<<present_p[3]<<endl;

 //=======paramters initialization=========


//============================================

	
	
}
