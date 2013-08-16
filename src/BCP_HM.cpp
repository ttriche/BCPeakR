#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include  <cassert>
#include	<string>
#include	<vector>
#include	<map>
#include	<cmath>
#include	<iomanip>
#include 	<algorithm>
#include	<stdlib.h>
#include	<stdio.h>
#include	<unistd.h>
#include  "TNT/tnt.h"
#include  "functions_IO.h"
#include	"functions_LU.h"
#include	"functions_Cholesky.h"
#include  "functions_Math.h"
#include  "cppoisson_HM.cpp"

#define   MINARGS 3

using namespace std;
using namespace TNT;

int num_win,num_data_frag,num_input_frag,num_seg,num_allseg;
double len_bf,len_aft,sum_bf,sum_aft,average_bf,average_aft,cutline;

const int L1=3000000;
const int L2=1500000;
const int N0=1;
const int N=26;
vector<int> Weight(L1);

vector< vector<int> > Plus_data(N);
vector< vector<int> > Minus_data(N);
vector< vector<int> > Plus_input(N);
vector< vector<int> > Minus_input(N);


double (*win_data)[4];
double (*ss)[6];
double (*rseg)[7];
int *data_frag;
int *input_frag;
double st[100001];

void printhelp();
int cmp (const void *a , const void *b );
void Add(int aa, int bb);
void trans2window(int r);
void seg(Matrix<double>,double);
int frag_count1(int);
int frag_count2(int);
double prob_pois(int,double);

int needhelp;
int next_option;
int frag_size = 200;
int win_size =200;
double p_value = 0.001;

char * ChIP_data = 0;
char * Input_data = 0;
char * Seg_Results = 0;

int main(int argc, char* argv[]) { // {{{

  // FIXME: test
  needhelp = 0;

	while ( (next_option = getopt(argc, argv, "1:2:f:w:p:3:h:")) != -1) {
    
    int this_option_optind = optind ? optind : 1;

    switch (next_option) {

			case '1': ChIP_data = optarg;  
                needhelp++;
                break;

      case '2': Input_data = optarg;
                needhelp++;
                break;

			case 'f': sscanf (optarg,"%d",&frag_size);
                break;

			case 'w': sscanf (optarg,"%d",&win_size);
                break;

			case 'p': sscanf (optarg,"%lf",&p_value); 
                break;

			case '3': Seg_Results = optarg;
                needhelp++;
                break;

			case 'h': sscanf (optarg,"%d",&needhelp);	
                break;

			default:  break;
		}
	}
  if(needhelp < 3) {
    printhelp();
  } else {	

    FILE *out1,*out3,*out4;
	
  	int *temp;
	  string filename,test_name,chrName, info;
    char chrSign,name4print[20];
    int m,chrNum,a,b,c1,pre,num,t1,t2,l1,l2,seg_len;
  	double num_temp, alpha, beta,thre,c,max,max1,lambda;

	filename = ChIP_data;

	num_win=0;


	ifstream fin(filename.c_str());
	if(!fin.is_open()){
		cout<<"The file "<<filename<<" doesn't exist."<<endl;
		exit(1);
	}

	while(!fin.eof())
	{
		fin>>chrName>>a>>b>>info>>c1>>chrSign;
		if(fin.fail()) break;

		test_name = chrName.substr(0,3);
		if (test_name!="chr") {
      cout<<"The ChIP-seq data is not in valid BED6 format."<<endl; 
      exit(1);
    }
	
    // FIXME: this will fail for non-human peak calling  
    // FIXME: code duplication, break out into a function
		chrName=chrName.substr(chrName.find('r')+1,chrName.size()-1);
		if(chrName=="x"||chrName=="X") chrNum=23;
		else if (chrName=="y" || chrName=="Y") chrNum=24;
		else if (chrName=="m" || chrName=="M") chrNum=25;
		else chrNum=atoi(chrName.c_str());

    // FIXME: code duplication, break out into a function
		if(chrNum<N) {
			if(chrSign=='+') {
        Plus_data[chrNum].push_back(a);
      } else if(chrSign=='-' && b>frag_size) {
        Minus_data[chrNum].push_back(b+1-frag_size);
      }
		}
	}
	fin.close();
	fin.clear();

	filename= Input_data;

	fin.open(filename.c_str());
	if(!fin.is_open()){
 		cout<<"The file "<<filename<<" doesn't exist."<<endl;
		exit(1);
	}

	while(!fin.eof())
	{
		fin>>chrName>>a>>b>>info>>c1>>chrSign;
		if(fin.fail()) break;

		test_name = chrName.substr(0,3);
		if (test_name!="chr") {
      // FIXME: code duplication; turn this into a function
      cout << "The control data is not in valid BED6 format." << endl; 
      exit(1);
    }
    // FIXME: should become test_bed6(test_name, 'control')
				
    // FIXME: code duplication; turn this into a function
		chrName=chrName.substr(chrName.find('r')+1,chrName.size()-1);
		if(chrName=="x"||chrName=="X") chrNum=23;
		else if (chrName=="y" || chrName=="Y") chrNum=24;
		else if (chrName=="m" || chrName=="M") chrNum=25;
		else chrNum=atoi(chrName.c_str());
    // FIXME: should become getChrNum(chrName)

    // FIXME: code duplication; turn this into a function
		if(chrNum<N){
			if(chrSign=='+') {
        Plus_input[chrNum].push_back(a);
      } else if(chrSign=='-' && b>frag_size) {
        Minus_input[chrNum].push_back(b+1-frag_size);
      }
		}
	}
	fin.close();


	out1 = fopen(Seg_Results,"w"); // FIXME: output BED instead of BCP (BED5)


	// FIXME: factor this out into a callable function, e.g.
	//
	// GRanges *BCP_eval( *ChIPRanges, *InputRanges, *options ) 
	//
	// presumably, this is where we want to start if using e.g. GRanges...
	for(int i=N0; i<N; i++)
	{

		/*The data part*/
		data_frag= new int[Plus_data[i].size()+Minus_data[i].size()+1];
		num=1;data_frag[0]=0;
		if(Plus_data[i].size()>0)
		{
			temp= new int[Plus_data[i].size()];

			for(unsigned int j=0; j<Plus_data[i].size(); j++)
				temp[j] = Plus_data[i][j];

			qsort(temp,Plus_data[i].size(),sizeof(int),cmp);
		
			pre=-1;
			for(unsigned int j=0; j<Plus_data[i].size(); j++)
			{
				if(temp[j]!=pre)
				{
					data_frag[num]=temp[j];
					Add(temp[j],temp[j]+frag_size-1);
					num++;
				}
				pre=temp[j];
			}
			delete []temp;
		}
		
		if(Minus_data[i].size()>0)
		{
			temp= new int[Minus_data[i].size()];
		
			for(unsigned int j=0; j<Minus_data[i].size(); j++)
				temp[j] = Minus_data[i][j];

			qsort(temp,Minus_data[i].size(),sizeof(int),cmp);
		
			pre=-1;
			for(unsigned int j=0; j<Minus_data[i].size(); j++)
			{
				if(temp[j]!=pre)
				{
					data_frag[num]=temp[j];
					Add(temp[j],temp[j]+frag_size-1);
					num++;
				}
				pre=temp[j];
			}
			delete []temp;
		}
		num_data_frag = num;// num_data should be num-1, just convenient for qsort;
		cout<<"number of data fragments:"<<num_data_frag<<endl;
		qsort(data_frag,num_data_frag,sizeof(int),cmp);

		//The input part
		input_frag= new int[Plus_input[i].size()+Minus_input[i].size()+1];
		num=1; input_frag[0] =0;
		if(Plus_input[i].size()>0)
		{
			temp= new int[Plus_input[i].size()];
		
			for(unsigned int j=0; j<Plus_input[i].size(); j++)
				temp[j] = Plus_input[i][j];

			qsort(temp,Plus_input[i].size(),sizeof(int),cmp);
		
			pre=-1;
			for(unsigned int j=0; j<Plus_input[i].size(); j++)
			{
				if(temp[j]!=pre)
				{
					input_frag[num]=temp[j];
					num++;
				}
				pre=temp[j];
			}
			delete []temp;
		}

		if(Minus_input[i].size()>0)
		{
			temp= new int[Minus_input[i].size()];
		
			for(unsigned int j=0; j<Minus_input[i].size(); j++)
				temp[j] = Minus_input[i][j];

			qsort(temp,Minus_input[i].size(),sizeof(int),cmp);
		
			pre=-1;
			for(unsigned int j=0; j<Minus_input[i].size(); j++)
			{
				if(temp[j]!=pre)
				{
					input_frag[num]=temp[j];
					num++;
				}
				pre=temp[j];
			}
			delete []temp;
		}
		num_input_frag = num;// num_input should be num-1, just convenient for qsort
		cout<<"number of input fragments:"<<num_input_frag<<endl;
		qsort(input_frag,num_input_frag,sizeof(int),cmp);

		cout<<"finished preprocessing data for chr"<<i<<endl;

		if(Plus_data[i].size()+Minus_data[i].size()>0)
		{	
			win_data= new double[L2][4];
			trans2window(i);			

			len_bf=sum_bf=0.0;
			for(m=2;m<=num_win;m++){
                		len_bf+= win_data[m][3];
                		sum_bf +=win_data[m][2]*win_data[m][3];
        		}
        		average_bf =sum_bf/len_bf;


			Matrix <double> obs(num_win+1,4);
			Matrix <double> data(num_win+1,5);	

			m=0;
			obs[m][0]=obs[m][1]=obs[m][2]=obs[m][3]=0.0; 
			for(m=1;m<=num_win;m++){
				obs[m][0] = win_data[m][0];
				obs[m][1] = win_data[m][1];
				obs[m][2] = floor(win_data[m][2]+0.5);
				obs[m][3] = win_data[m][3];
			}

			len_aft=sum_aft=0.0;
			for(m=2;m<=num_win;m++){
				len_aft+= obs[m][3];
				sum_aft +=obs[m][2]*obs[m][3];
			}
			average_aft =sum_aft/len_aft;

      // FIXME: terrifying
			for(m=1;m<num_win;m++){
				if((obs[m][2]==0)&&(obs[m][3]==1)&&((obs[m-1][2]!=0)||(obs[m+1][2]!=0)))
          obs[m][3]=0;
      }

			num_temp = floor(log10((double)num_win)+0.5);
			double p = 1.0/pow(10.0,num_temp);
			if(average_aft<0.1&&p<0.0001) p=p*100;
			if(average_aft>=0.1&&average_aft<=0.5&&p<0.0001) p=p*10;
     
			if (win_size>=200) {
        t1=10;
        t2=5;
      }
			if (win_size<200 && win_size>=100){
        t1=7;
        t2=4;
      }
			if (win_size<100){
        t1=5;
        t2=3;
      }
			if (t1>=num_win) {
        t1=5;
        t2=3;
      }
			if (t1>=num_win) {
        cout<<"Too little data for the model."<<endl; 
        exit(1);
      }

			cppoisson  tmp(obs, p, 1.0, 1.0, t1, t2);
			tmp.BcmixSmooth(); // horrifying
			
			m=0;
			data[m][0]=data[m][1]=data[m][2]=data[m][3]=data[m][4]=0.0;
        		for(m=1;m<=num_win;m++){
                		data[m][0] = obs[m][0];
                		data[m][1] = obs[m][1];
                		data[m][2] = obs[m][2];
                		data[m][3] = obs[m][3];
                		data[m][4] = tmp.estPara[m];
        		}
			cout<<"finished calculating posterior means of chr"<<i<<endl;
			           
			thre = 0.9;
			//calculating the factorial first
			st[0]=0;
			for (m=1;m<=100000;m++)
				st[m]= st[m-1]+log(m);
			
			m=0;
			while (prob_pois(m,average_aft)<thre) m++;

			if (m==0) cutline = 1.0;
			if (m>3) cutline = m*1.0-0.5;
			if (m>=1&&m<=3) cutline = m*1.0;


			ss = new double [num_win][6];

			seg(data,cutline);
         
			rseg = new double[num_seg+1][7];

			int m1=0;
			for (m=1;m<=num_allseg;m++){
				if(ss[m][3]>=cutline){                             
					m1++;
					rseg[m1][0]=ss[m][0];
					rseg[m1][1]=ss[m][1];
					rseg[m1][2]=ss[m][2];
					rseg[m1][3]=ss[m][3];
					rseg[m1][4]=rseg[m1][5]=0;
					rseg[m1][6]=1.0;
				}
			}
			
			c=(double)num_data_frag/(double)num_input_frag; // horrifying

			max=(num_input_frag-1)/(double)(input_frag[num_input_frag-1]+frag_size-1-input_frag[1]);

		
			for (m=1;m<=num_seg;m++){
				seg_len = (int)rseg[m][2];
				max1 = max*seg_len;
			
				l1 = frag_count1(m);
				rseg[m][4] = (double)l1;

				l2 = frag_count2(m);
				lambda=l2>=max1?l2:max1;

				rseg[m][5] = c*lambda;
				rseg[m][6] = fabs(1.0-prob_pois(l1,rseg[m][5]));
			}

      // FIXME: should be a function
			if (i>=1&&i<=22) sprintf(name4print,"chr%d",i);
			if (i==23) sprintf(name4print,"chrX");
			if (i==24) sprintf(name4print,"chrY");
			if (i==25) sprintf(name4print,"chrM");

      for (m=1;m<=num_seg;m++) {
        if((rseg[m][6]<=p_value)&&(rseg[m][4]>rseg[m][5])) {

          // FIXED: switch to BED6
  				fprintf(out1, "%s\t%d\t%d\t%c\t%lf\t%c\n",
                        name4print,           // chrom
                        (int)(rseg[m][0]-1),  // start
                        (int)rseg[m][1],      // end
                        '.',                  // name
                        rseg[m][3],           // score
                        '*'                   // strand
                        ); // BED6 output!
        }
      }

			delete [] rseg;
			delete [] ss;
			delete [] win_data;

			cout<<"Finished choosing significant segments of chr"<<i<<endl;

		} else {
      cout<<"chr"<<i<<"is empty."<<endl;
    }
	
		for(m=0; m< L1; m++) Weight[m]=0;     
		delete [] data_frag;
		delete [] input_frag;
	}
	fclose(out1);
  }
	return 0;
} // }}}

void printhelp() { // {{{
    cout<<""<<endl;
    cout<<"BCP_HM: a Bayesian change-point based peak caller for histone marks."
        <<endl;
    cout<<""<<endl;
    cout<<"Usage:"<<endl;
    cout<<""<<endl;
    cout<<"\tBCP_HM -1 chipseq.bed -2 input.bed -3 results.bed"<<endl;
    cout<<""<<endl;
    cout<<"The following are required:"<<endl;
    cout<<""<<endl;
    cout<<"\t-1\tThe ChIP-seq data set you want to call peaks from."<<endl;
    cout<<"\t-2\tThe control/input data set you want to use."<<endl;
    cout<<"\t-3\tThe filename for peak calls to be output."<<endl;
    cout<<""<<endl;
    cout<<"The following are optional:"<<endl;
    cout<<""<<endl;
	  cout<<"\t-f\tThe fragment size to which we extend reads in preprocessing."
        <<endl;
    cout<<"\t\t(Default is 200bp)"<<endl;
  	cout<<"\t-w\tThe window size we apply to adjacent windows in preprocessing."
        <<endl;
    cout<<"\t\t(Default is 200bp)"<<endl;
    cout<<"\t-p\tThe p_value cutoff to avoid false positives based on input."
        <<endl;
    cout<<"\t\t(Range: 1e-2 to 1e-6, default is 1e-3)"<<endl;
    cout<<""<<endl;
} // }}}

void Add(int aa, int bb) { // {{{
	aa++;
	int head_index=aa/win_size;
	int tail_index=bb/win_size;

	
	if(head_index==tail_index)	Weight[head_index] += bb-aa+1;
	else //if(head_index<tail_index)
	{
		
		Weight[head_index] += win_size - aa%win_size;

		for(int k=1; k<tail_index-head_index; k++)
			Weight[head_index+k] += win_size;

		Weight[tail_index] += bb%win_size+1;

	}
} // }}}
			
void trans2window(int r) { // {{{
	int pre_index=-1;
	int i1=0;

	for(unsigned int j=0;j<L1;j++)
	{
		if(Weight[j]>0)
		{
			int n1=j-pre_index-1;
			if(n1>0)
			{
				i1++;
				win_data[i1][0]=(pre_index+1)*win_size;
				win_data[i1][1]=win_data[i1][0]+n1*win_size-1;
				win_data[i1][2]=0;
				win_data[i1][3]=n1;
			}

			i1++;
			win_data[i1][0]=j*win_size;
			win_data[i1][1]=win_data[i1][0]+win_size-1;
			win_data[i1][2]=Weight[j]*1.0/win_size;
			win_data[i1][3]=1;

			pre_index=j;
		}
	}
	num_win=i1;
	cout<<"number of window of chr"<<r<<"\t"<<num_win<<endl;

} // }}}

int cmp ( const void *a , const void *b ) { // {{{
return *(int *)a - *(int *)b;
} // }}}

double prob_pois(int l1, double lambda){ // {{{
        int i2;
        double prob,logprob,sum_temp;
        prob = exp(-1*lambda);
        for (i2=1;i2<=l1;i2++){
                logprob=0.0;
                sum_temp = st[i2];
                logprob = (-lambda)+i2*log(lambda)-sum_temp;
                prob += exp(logprob);
        }
        return prob;
} // }}}

void  seg(Matrix<double> input,double thre){ // {{{
        int i,j,k,j1;
        double temp_sum;
        j=1;i=1;num_seg=0;

        while (1){
                temp_sum = 0.0;
                if (j> num_win) break;
                if (input[j][4] <thre){
                        j1=j;
                        while ((j1<=num_win)&&(input[j1][4] < thre)) j1++;
                        ss[i][0] = input[j][0];
                        ss[i][1] = input[j1-1][1];
                        ss[i][2] = ss[i][1] -ss[i][0] + 1.0;
                        for (k=j; k<j1; k++) 
                                temp_sum += (input[k][1]-input[k][0]+1.0)*input[k][4];
                        ss[i][3] = temp_sum/ss[i][2];
                        ss[i][4] = j;
                        ss[i][5] = j1-1;
                        i++;
                        j = j1;
                }
                else {
			num_seg++;
                        j1=j;
                        while ((j1<=num_win)&&(input[j1][4] >= thre)) j1++;
                        ss[i][0] = input[j][0];
                        ss[i][1] = input[j1-1][1];
                        ss[i][2] = ss[i][1] -ss[i][0] + 1.0;
                        for (k=j; k<j1; k++) 
                                temp_sum += (input[k][1]-input[k][0]+1.0)*input[k][4];
                        ss[i][3] = temp_sum/ss[i][2];
                        ss[i][4] = j;
                        ss[i][5] = j1-1;
                        i++;
                        j = j1;

                }
        }
        num_allseg = i-1;
	//cout<<"number of all segs:"<<num_allseg<<endl;
	cout<<"number of candidate segments:"<<num_seg<<endl;
} // }}}

int frag_count1(int rec){ // {{{

        int start_ind,end_ind,comp_ind;
        int ind1,ind2;
        int count;

        start_ind =1;
        end_ind = num_data_frag-1;
        while (start_ind<end_ind-1){
                comp_ind=start_ind+(end_ind-start_ind)/2;
                if(data_frag[comp_ind]==rseg[rec][1]){
			start_ind = comp_ind;
			break;
		}
		else{
			if(data_frag[comp_ind]<rseg[rec][1])
				start_ind=comp_ind;
			else
				end_ind=comp_ind;
		}

	}
	ind1 = start_ind;


	start_ind =1;
	end_ind = num_data_frag-1;
	while(start_ind<end_ind-1){
		comp_ind=start_ind+(end_ind-start_ind)/2;
		if((data_frag[comp_ind]+frag_size-1) ==rseg[rec][0]){
			end_ind = comp_ind;
			break;
		}
		else{
			if((data_frag[comp_ind]+frag_size-1)>rseg[rec][0])
				end_ind=comp_ind;
			else start_ind = comp_ind;
		}
	}
	ind2=end_ind;

	count = ind1-ind2+1;
	return count;
} // }}}

int frag_count2(int rec){ // {{{
        int start_ind,end_ind,comp_ind;
        int ind1,ind2;
        int count;

        start_ind =1;
        end_ind = num_input_frag-1;
        while (start_ind<end_ind-1){
                comp_ind=start_ind+(end_ind-start_ind)/2;
                if(input_frag[comp_ind]==rseg[rec][1]){
                        start_ind = comp_ind;
                        break;
                }
                else{
                        if(input_frag[comp_ind]<rseg[rec][1])
                                start_ind=comp_ind;
                        else
                                end_ind=comp_ind;
                }
                
        }
        ind1 = start_ind;


        start_ind =1;
        end_ind = num_input_frag-1;
        while(start_ind<end_ind-1){
                comp_ind=start_ind+(end_ind-start_ind)/2;
                if((input_frag[comp_ind]+frag_size-1) ==rseg[rec][0]){
                        end_ind = comp_ind;
                        break;
                }
                else{
                        if((input_frag[comp_ind]+frag_size-1)>rseg[rec][0])
                                end_ind=comp_ind;
                        else start_ind = comp_ind;
                }
        }       
        ind2=end_ind;

        count = ind1-ind2+1;
        return count;
} // }}}
