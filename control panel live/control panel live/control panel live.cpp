// control panel live.cpp : Defines the entry point for the console application.
//



#include "stdafx.h"
#include<stdio.h>
#include <iostream>
#include<conio.h>
#include<fstream>
#include<cmath>
#include<math.h>
#include<string>
#include<Windows.h>
using namespace std;

#define N 5
#define M 32
#define T 85

int accuracy=0;							//TO CALCULATE ACCURACY
int xander=0;							//TO STORE VALUES IN PROBABILITY_SEQ ARRAY
long double initial_samples[150000];    //TO STORE RECORDING DATA
long double noise[150000];              //TO STORE NOISE DATA
long double stable[10400];				//TO STORE 7040 OVERLAPPED STABLE SPEECH SAMPLES
int count_samples=0;					//TO STORE COUNT OF INITIAL_SAMPLES ARRAY
int index;								//TO STORE INDEX TO FIND STABLE SPEECH SAMPLES
long double final_ci[85][12];			//TO STORE CI VALUES OF 85 FRAMES OF ONE RECORDING
int obs_seq[85];						//TO STORE OBSERVATION SEQUENCE
long double universe_codebook[32][12];  //TO STORE REFERENCE CODEBOOK
long double raise_sine[12]={2.552914271,4,5.242640687,6.196152423,6.795554958,7,6.795554958,6.196152423,5.242640687,4,2.552914271,1};   //TO STORE RAISED SINE WINDOW
long double tok_weights[12]={1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0};       //TO STORE TOKURA WEIGHTS
long double hamming_window[320];        //TO STORE HAMMING WINDOW 

long double A[N][N];              //TO STORE A MATRIX
long double B[N][M];			  //TO STORE B MATRIX
long double pi[N];				  //TO STORE PI MATRIX
long double alpha[T][N];          //TO STORE ALPHA MATRIX
long double beta[T][N];           //TO STORE BETA MATRIX
long double gamma[T][N];          //TO STORE GAMMA MATRIX
long double delta[T][N];          //TO STORE DELTA MATRIX
int si[T][N];                     //TO STORE SI MATRIX
long double pstar;                //TO STORE PSTAR VALUE
int qstar;                        //TO STORE QSTAR MATRIX
int path[T];                      //TO STORE STATE SEQUENCE
long double xi[T][N][N];          //TO STORE XI MATRIX
long double pi_bar[N];            //TO STORE UPDATED PI MATRIX
long double A_bar[N][N];          //TO STORE UPDATED A MATRIX
long double B_bar[N][M];          //TO STORE UPDATED B MATRIX
long double temp_A[N][N];         //TO CALCULATE AVERAGE A MATRIX
long double temp_B[N][M];         //TO CALCULATE AVERAGE B MATRIX
long double pass_A[N][N];		  //TO PASS AVERAGE A MATRIX TO NEXT ITERATION OF HMM
long double pass_B[N][M];         //TO PASS AVERAGE B MATRIX TO NEXT ITERATION OF HMM
long double prob_seq[10];         //TO STORE PROBABLITY OF OBS_SEQ WITH 10 DIGITS

void input_universal_codebook();               //TO INPUT REFERENCE CODEBOOK
void input(int,int);                           //TO INPUT 200 TRAINING RECORDINGS
void input_test(int);                          //TO INPUT 100 TESTING RECORDINGS
void do_dc_shift();                            //TO CALCULATE DC SHIFT
void Normalise();                              //TO NORMALISE DATA
void select_frames();                          //TO SELECT 85 OVERLAPPED STABLE FRAMES
void durbin(int,long double[],long double[]);  //TO FIND AI
void calc_ci(long double[],long double[]);     //TO CALCULATE CI
void find_observation_sequence();              //TO FIND OBSERVATION SEQUENCE
void tokura();                                 //TO CALCULATE TOKURA DISTANCE
void input_hamming_window();                   //TO INPUT HAMMING WINDOW
void input_test(int,int);                      //TO INPUT 100 TESTING RECORDINGS

void input_model();
void input_A();				//TO INPUT MATRIX A
void input_B();				//TO INPUT MATRIX B
void input_pi();			//TO INPUT PI
void forward(int[]);		//TO RUN FORWARD PROCEDURE
void backward(int[]);	    //TO RUN BACKWARD PROCEDURE
void move_to_text(int);     //TO MOVE ALPHA,BETA,GAMMA,DELTA,SI TO TEXT FILE
void calc_gamma();          //TO CALCULATE GAMMA
void viterbi(int[]);        //TO CALCULATE DELTA AND SI
void reestimation(int[]);   //TO CALCULATE XI,PI_BAR,A_BAR,B_BAR
void initialize_temp_A_B(); //TO MAKE TEMP_A AND TEMP_B ARRAY TO 0
void input_model2();        //TO INPUT RECENT A AND B MATRIX
void testing();             //TO TEST DATA
void read_digit_model();    //TO READ FINAL MODELS OF EACH DIGIT

void live_train();
void live_input(int);
void live_move_to_text();

int main()
{	
	pi[0]=1;
	for(int j=1;j<N;j++)
	{
		pi[j]=0;
	}
	input_universal_codebook();
	input_hamming_window();
	
	//TRAINING STARTS HERE
	/*
	for(int j=0;j<10;j++)
	{
		//cout<<"DIGIT "<<j<<" TRAINING : "<<endl;
		for(int count=0;count<3;count++)
		{
			//cout<<endl<<"DIGIT "<<j<<" "<<"RUN NO: "<<count+1<<endl;
			initialize_temp_A_B();
			//FOR 20 UTTERANCE OF SAME DIGIT
			for(int i=0;i<20;i++)
			{
				input(i,j);
				do_dc_shift();
				Normalise();
				select_frames();       
				find_observation_sequence();

				if(count==0)
				{
					 input_model();  //TO INPUT THE MODEL AND OBSERVATION SEQUENCES
				}
				else
				{
					input_model2();
				}
		      
				int turn=0;
				while(turn!=20)
				{
					//cout<<"*******************************************************************"<<endl;
					//cout<<"Iteration Number : "<<turn+1<<endl;
					forward(obs_seq);   //TO RUN FORWARD PROCEDURE
					//forward(obs_seq);  //TO RUN FORWARD PROCEDURE
					backward(obs_seq); //TO RUN BACKWARD PROCEDURE
					calc_gamma();    //TO CALCULATE GAMMA
					viterbi(obs_seq);  //TO CALCULATE DELTA AND SI
					reestimation(obs_seq);
					//cout<<endl;
					turn++;
				}

				//KEEP ON ADDING TO TEMP_A AND TEMP_B
				for(int r=0;r<N;r++)
				{
					for(int t=0;t<N;t++)
					{
						temp_A[r][t]=temp_A[r][t]+A[r][t];
						//cout<<temp_A[r][t]<<"  ";
					}
					//cout<<endl;
				}

				for(int r=0;r<N;r++)
				{
					for(int t=0;t<M;t++)
					{
						temp_B[r][t]=temp_B[r][t]+B[r][t];
						//cout<<temp_B[r][t]<<"  ";
					}
					//cout<<endl;
				}
			}

			//DIVIDE TEMP_A AND TEMP_B BY 20
			for(int r=0;r<N;r++)
				{
					for(int t=0;t<N;t++)
					{
						temp_A[r][t]=temp_A[r][t]/20;
						pass_A[r][t]=temp_A[r][t];
						//cout<<temp_A[r][t]<<"  ";
					}
					//cout<<endl;
				}

			for(int r=0;r<N;r++)
				{
					for(int t=0;t<M;t++)
					{
						temp_B[r][t]=temp_B[r][t]/20;
						pass_B[r][t]=temp_B[r][t];
						//cout<<temp_B[r][t]<<"  ";
					}
					//cout<<endl;
				}
		}

		//cout<<endl<<endl;
	
	//STORE FINAL MODELS
	move_to_text(j);  //TO MOVE ALPHA,BETA,GAMMA,DELTA,SI,XI TO TEXT FILE
	}
	*/

	//TRAINING ENDS HERE
	//cout<<"TRAINING DONE"<<endl;

	//testing();

	live_train();

	getch();
	return 0;
}


//READ GIVEN UNIVERSAL CODEBOOK
void input_universal_codebook()
{
	long double data;
	int i=0,j=0,count=0;
	ifstream obj("codebook.txt");
	while(obj>>data)
	{
		if(j==12)
		{
			j=0;
			i++;
		}
		universe_codebook[i][j]=data;
		j++;
		count++;
	}
}

//TO INPUT HAMMING WINDOW
void input_hamming_window()
{
	ifstream obj("windowFunction.txt");
	long double data;
	int z=0;
	while(obj>>data)
	{
		hamming_window[z]=data;
		z++;
	}
}

//TO MAKE TEMP_A & TEMP_B ARRAY TO 0
void initialize_temp_A_B()
{
	for(int r=0;r<N;r++)
		{
			for(int t=0;t<N;t++)
			{
				temp_A[r][t]=0;
			}
		}

	for(int r=0;r<N;r++)
		{
			for(int t=0;t<M;t++)
			{
				temp_B[r][t]=0;
			}
		}
}

//INPUTS THE TRAINING FILES
void input(int i,int j)
{
	long long int h=i+1;
	long long int f=j;
	long double data1;
	int z=0;
	count_samples=0;
	ifstream obj;

	string filename;
		filename="194101013_"+to_string(f)+"_"+to_string(h)+".txt";
	
	obj.open(filename);
	
	while(obj>>data1)
	{
		initial_samples[z]=data1;
		//cout<<initial_samples[z]<<endl;
		z++;
	}
	count_samples=z;
	//cout<<"number of samples : "<<count_samples<<endl;
}


//CALCULATES DC_SHIFT
void do_dc_shift()
{
	long double data2;
	long double sum=0;
	long double dc_shift=0;
	int z=0;

	ifstream obj("Silence.txt");
	while(obj>>data2)
	{
		noise[z]=data2;
		sum=sum+data2;
		z++;
	}
	dc_shift=sum/z;
	//cout<<"dc_shift is : "<<dc_shift<<endl;

	//SUBSTRACT THE DC SHIFT VALUE FROM EACH SAMPLE 
	for(int i=0;i<count_samples;i++)
	{
		initial_samples[i]=initial_samples[i]-dc_shift;
	}		
}


//NORMALISE DATA IN A RANGE OF [-5000,5000]
void Normalise()
{
	long double max=0;
	long double min=0;
	long double normalise;
	int index_max=0,index_min=0;

	//FIND MAX AND MIN VALUE
	for(int i=0;i<count_samples;i++)
	{
		if(initial_samples[i]>max)
		{
			max=initial_samples[i];
			index_max=i;
		}
		if(initial_samples[i]<min)
		{
			min=initial_samples[i];
			index_min=i;
		}
	}
	//cout<<"max: "<<max<<endl;
	//cout<<"min: "<<min<<endl;

	//DECIDE WHICH ABSOLUTE VALUE IS MORE. USE THAT FOR NORMALISATION
	if(fabs(max)>=fabs(min))
	{
		normalise=fabs(max);
		index=index_max;
	}
	else
	{
		normalise=fabs(min);
		index=index_min;
	}
	//cout<<"parameter used for normalisation is : "<<normalise<<endl;
	//cout<<"index "<<index;

	//NORMALISE USING THE FOUNDED VALUE
	for(int i=0;i<count_samples;i++)
	{
		initial_samples[i]=(initial_samples[i]/normalise)*5000;
	}
}



//SELECT 85 FRAMES AND CALCULATES ITS CEPSTRAL COEFFICIENTS
void select_frames()
{
	int z=0,k=0,f=0;
	
	//SELECTING 85 OVERLAPPED FRAMES
	for(int i=(index-5200);i<=(index+5199);i++)
	{
		stable[z]=initial_samples[i];
		//cout<<stable[z]<<endl;
		z++;
	}
	//cout<<"Number of samples in stable frames : "<<z<<endl;

	//FOR 85 FRAMES CALCULATE Ri Ai Ci
	for(int i=0;i<85;i++)
	{
		long double samples[320];
		int p=12;
		long double R[13];
		long double a[13];
		long double c[13];
		long double val=0;

		//cout<<"frame : "<<i<<endl;
		//READING 320 SAMPLES INTO AN ARRAY 
		for(int j=0;j<320;j++)
		{
			samples[j]=stable[j+k]*hamming_window[j];
			//cout<<samples[j]<<endl;
		}
		//cout<<endl<<endl;
		k=k+120;
		

		//CALCULATING Ri
		int m=0,n=0;
		//cout<<"values of Ri are "<<endl;
		for(n=0;n<=12;n++)
		{
			for(m=0;m<320-n;m++)
			{
				val=val+samples[m]*samples[m+n];
			}
			R[n]=val;
			//cout<<"R"<<n<<"  "<<R[n]<<endl;
			val=0;
		}

		//FUNCTION CALL TO durbin...
		durbin(p,R,a);

		//FUNCTION CALL TO calc_ci...
		calc_ci(a,c);
		
		//STORING ALL Ci VALUES IN A 2D ARRAY
		for(int x=0;x<12;x++)
		{
			final_ci[f][x]=c[x+1]*raise_sine[x];
			//cout<<"ci["<<f<<","<<x<<"]   "<<final_ci[f][x]<<endl;
		}
		f++;	
	}
}


//THIS FUNCTION WILL GIVE Ai VALUES AS OUTPUT
void durbin(int p,long double R[],long double a[])
{
    long double E[13];
    long double K[13];
    long double alpha[13][13];
    long double internal_val;

	//STEP 1: CALCULATE E[0]
    E[0]=R[0];
    int j,i;


	//STEP 2: CALCULATE K[i]
    for(i=1;i<=p;i++)
    {
        if(i==1)
        {
            internal_val=0;
        }
        for(j=1;j<=i-1;j++)
        {
            internal_val=internal_val+alpha[j][i-1]*R[i-j];
        }
        K[i]=(R[i]-internal_val)/E[i-1];
		internal_val=0;
		

		//STEP 3: CALCULATE ALPHA[i][i]; 
        alpha[i][i]=K[i];


		//CALCULATE ALPHA[j][i]
        for(j=1;j<=i-1;j++)
        {
            alpha[j][i]=alpha[j][i-1]-(K[i]*alpha[i-j][i-1]);
        }

		//STEP 5: CALCULATE E[i]
        E[i]=(1-K[i]*K[i])*E[i-1];
	}

	//PRINT VALUES OF Ai
    //cout<<"values of ai are"<<endl;
    for(i=1;i<=p;i++)
    {
        a[i]=alpha[i][p];
        //cout<<"a"<<i<<"   "<<a[i]<<endl;
    }
}


//THIS FUNCTION WILL GIVE Ci VALUES AS OUTPUT
void calc_ci(long double a[],long double c[])
{
	int i=0,j=0;
	long double internal_val=0;
	long double k;
	for(i=1;i<=12;i++)
	{
		if(i==1)
		{
			internal_val=0;
		}
		for(j=1;j<=i-1;j++)
		{
			k=long double(j)/long double(i);
			internal_val=internal_val+(k*c[j]*a[i-j]);
		}
		c[i]=a[i]+internal_val;
		internal_val=0;
	}
}

//FIND OBSERVATION SEQUENCE
void find_observation_sequence()
{
	for(int i=0;i<85;i++)
	{
		long double distance[32];
		for(int j=0;j<32;j++)
		{
			distance[j]=0;
		}

		for(int j=0;j<32;j++)
		{
			long double temp=0;
			for(int k=0;k<12;k++)
			{
				temp=(final_ci[i][k]-universe_codebook[j][k])*(final_ci[i][k]-universe_codebook[j][k]);
				distance[j]=distance[j]+(tok_weights[k]*temp);
			}
		}

		long double min=INT_MAX;
		int indx=0;
		for(int y=0;y<32;y++)
		{
			if(distance[y]<min)
			{
				min=distance[y];
				indx=y;
			}
		}
		obs_seq[i]=indx+1;
		//cout<<obs_seq[i]<<"  ";

		/*for(int j=0;j<32;j++)
		{
			//cout<<"distance"<<j<<" "<<distance[j]<<endl;
		}
		cout<<endl;
		*/
	}
}


//TO INPUT THE MODEL AND OBSERVATION SEQUENCES 
void input_model()
{
	input_A();   //TO INPUT MATRIX A
	input_B();   //TO INPUT MATRIX B
	input_pi();  //TO INPUT PI
}

//TO UPDATE VALUES OF INITIAL MODEL A AND B WHILE DOING AVERAGING
void input_model2()
{
	for(int r=0;r<N;r++)
		{
			for(int t=0;t<N;t++)
			{
				A[r][t]=pass_A[r][t];
			}
		}

	for(int r=0;r<N;r++)
		{
			for(int t=0;t<M;t++)
			{
				B[r][t]=pass_B[r][t];
			}
		}
}

//TO INPUT MATRIX A
void input_A()
{
	ifstream obj;
		obj.open("A_matrix.txt");
	
	long double data;
	
	for(int j=0;j<N;j++)
	{
		for(int k=0;k<N;k++)
		{
			obj>>A[j][k];
	
			//cout<<A[j][k]<<"  ";
		}
		//cout<<endl;
	}
}

//TO INPUT MATRIX B
void input_B()
{
	ifstream obj("B_matrix.txt");
	long double arr[10000];
	long double data;
	int i=0;
	while(obj>>data)
	{
		arr[i]=data;
		i++;
	}
	i=0;
	for(int j=0;j<N;j++)
	{
		for(int k=0;k<M;k++)
		{
			B[j][k]=arr[i];
			i++;
			//cout<<B[j][k]<<"  ";
		}
		//cout<<endl;
	}
}

//TO INPUT PI
void input_pi()
{
	ifstream obj("pi.txt");
	long double data;
	int i=0;
	while(obj>>data)
	{
		pi[i]=data;
		//cout<<pi[i]<<"  ";
		i++;
	}
}



//TO RUN FORWARD PROCEDURE
void forward(int obs[])
{
	//INITIALIZATION
	for(int i=0;i<N;i++)
	{
		alpha[0][i]=pi[i]*B[i][obs[0]-1];
		//cout<<alpha[0][i]<<endl;
	}

	//INDUCTION
	for(int t=1;t<T;t++)
	{
		for(int j=0;j<N;j++)
		{
			long double sum=0;
			for(int i=0;i<N;i++)
			{
				sum=sum+(alpha[t-1][i]*A[i][j]);
				//cout<<"sum"<<sum<<endl;
			}
			alpha[t][j]=sum*B[j][obs[t]-1];
			//cout<<"alpha "<<t<<" "<<j<<" : "<<alpha[t][j]<<endl;
		}
	}

	//TERMINATION
	long double sum=0;
	for(int i=0;i<N;i++)
	{
		//cout<<alpha[T-1][i]<<"  ";
		sum=sum+alpha[T-1][i];
	}
	prob_seq[xander]=sum;
	//cout<<"probability of observation sequence : "<<sum<<endl;
}


//TO RUN BACKWARD PROCEDURE
void backward(int obs[])
{
	//INITIALIZATION
	for(int i=0;i<N;i++)
	{
		beta[T-1][i]=1;
	}

	//INDUCTION
	for(int t=T-2;t>=0;t--)
	{
		for(int i=0;i<N;i++)
		{
			long double sum=0;
			for(int j=0;j<N;j++)
			{
				sum=sum+(A[i][j]*B[j][obs[t+1]-1]*beta[t+1][j]);
			}
			beta[t][i]=sum;
			//cout<<beta[t][i]<<"  ";
		}
		//cout<<endl;
	}
}


//TO MOVE FINAL MODELS OF EACH DIGIT TO TEXT FILE
void move_to_text(int j)
{
	//TO MOVE NEW MODEL TO TEXT FILE
	long long int h=j;
	ofstream obj;
	string filename;
		filename="final_A_"+to_string(h)+".txt";
	obj.open(filename);

	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			obj<<temp_A[i][j]<<"    ";
		}
		obj<<endl;
	}

	ofstream obj2;
	string filename2;
		filename2="final_B_"+to_string(h)+".txt";
	obj2.open(filename2);

	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			obj2<<temp_B[i][j]<<"    ";
		}
		obj2<<endl;
	}
}


//TO CALCULATE GAMMA
void calc_gamma()
{
	for(int t=0;t<T;t++)
	{
		long double temp=0;
		for(int i=0;i<N;i++)
		{
			temp=temp+alpha[t][i]*beta[t][i];
		}
		for(int i=0;i<N;i++)
		{
			gamma[t][i]=(alpha[t][i]*beta[t][i])/temp;
			//cout<<gamma[t][i]<<"  ";
		}
		//cout<<endl;
	}
}


//TO CALCULATE DELTA AND SI
void viterbi(int obs[])
{
	//INITIALIZATION
	for(int i=0;i<N;i++)
	{
		delta[0][i]=pi[i]*B[i][obs[0]-1];
		si[0][i]=0;
	}

	//RECURSION
	for(int t=1;t<T;t++)
	{
		for(int j=0;j<N;j++)
		{
			long double max=0;
			int index=0;
			long double temp;
			for(int i=0;i<N;i++)
			{
				temp=delta[t-1][i]*A[i][j];
				if(temp>max)
				{
					max=temp;
					index=i;
				}
			}
			delta[t][j]=max*B[j][obs[t]-1];
			si[t][j]=index;
		}
	}

	//TERMINATION
	long double max=0;
	int index=0;
	long double temp;
	for(int i=0;i<N;i++)
	{
		temp=delta[T-1][i];
		if(temp>max)
		{
			max=temp;
			index=i;
		}
	}

	//CALCULATE PSTAR,QSTAR
	pstar=max;
	qstar=index;
	//cout<<"pstar : "<<pstar<<endl;
	//printf("pstar :%.200Lf",pstar);
	//cout<<"qstar : "<<qstar<<endl;

	//TO FIND STATE SEQUENCE
	path[T-1]=qstar;
	for(int t=T-2;t>=0;t--)
	{
		path[t]=si[t+1][path[t+1]];
	}

	/*
	//PRINT STATE SEQUENCE
	cout<<"STATE SEQUENCE"<<endl;
	for(int t=0;t<T;t++)
	{
		cout<<path[t]<<" ";
	}
	cout<<endl;
	*/
}

//CALCULATE UPDATED A AND B MATRIX
void reestimation(int obs[])
{
	//CALCULATE XI MATRIX
	for(int t=0;t<T-1;t++)
	{
		long double temp=0;
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				temp=temp+alpha[t][i]*A[i][j]*B[j][obs[t+1]-1]*beta[t+1][j];
			}
		}
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<N;j++)
			{
				xi[t][i][j]=(alpha[t][i]*A[i][j]*B[j][obs[t+1]-1]*beta[t+1][j])/temp;
				//cout<<xi[t][i][j]<<" ";
			}
			//cout<<endl;
		}
	}

	//MAKE LAST ROW OF XI 0
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			xi[T][i][j]=0;
		}
	}

	//TO CALCULATE PI BAR
	for(int i=0;i<N;i++)
	{
		pi_bar[i]=gamma[0][i];
		//cout<<pi_bar[i]<<"  ";
	}

	//COPY PI_BAR TO PI
	for(int i=0;i<N;i++)
	{
		pi[i]=pi_bar[i];
		//cout<<pi[i]<<"  ";
	}

	//TO CALCULATE A_BAR
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			long double temp1=0;
			long double temp2=0;
			for(int t=0;t<T-1;t++)
			{
				temp1=temp1+xi[t][i][j];
				temp2=temp2+gamma[t][i];
			}
			A_bar[i][j]=temp1/temp2;
			//cout<<A_bar[i][j]<<"  ";
		}
		//cout<<endl;
	}
	
	//COPY A_BAR TO A
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			A[i][j]=A_bar[i][j];
			//cout<<A[i][j]<<"  ";
		}
		//cout<<endl;
	}

	//TO CALCULATE B_BAR
	long double threshold;
	long double val1=10;
	long double val2=-30;
	threshold=pow(val1,val2);

	for(int j=0;j<N;j++)
	{
		int count=0;
		long double max=0;
		int index=0;
		for(int k=0;k<M;k++)
		{
			long double temp1=0;
			long double temp2=0;
			for(int t=0;t<T;t++)
			{
				if((obs[t]-1)==k)
				{
					temp1=temp1+gamma[t][j];
				}
				temp2=temp2+gamma[t][j];
			}
			B_bar[j][k]=temp1/temp2;

			//cout<<B_bar[j][k]<<"  ";

			if(B_bar[j][k]>max)
			{
				max=B_bar[j][k];
				index=k;
			}

			if(B_bar[j][k]<threshold)
			{
				B_bar[j][k]=threshold;
				count++;
			}
		}
		//cout<<endl;
		B_bar[j][index]=max-count*threshold;
	}

	//COPY B_BAR TO B
	for(int j=0;j<N;j++)
	{
		for(int k=0;k<M;k++)
		{
			B[j][k]=B_bar[j][k];
			//cout<<B[j][k]<<"  ";
		}
		//cout<<endl;
	}
}

//TO TEST FOR LIVE RECORDING
void testing()
{
		xander=0;
		
		long double data1;
		int z=0;
		count_samples=0;
		//RECORD YOUR VOWEL
		system("Recording_Module.exe 3 input_file.wav input_file.txt");
		ifstream obj1("input_file.txt");
		while(obj1>>data1)
		{
			initial_samples[z]=data1;
			//cout<<initial_samples[z]<<endl;
			z++;
		}
		count_samples=z;
		cout<<"number of samples for testing : "<<count_samples<<endl;

		do_dc_shift();
		Normalise();
		select_frames();       
		find_observation_sequence();
		
			//TO CALCULATE PROBABILITY WITH EACH DIGIT
		//1
			for(int u=0;u<8;u++)
			{
				int i=0,j=0;
				long long int s=u;
				ifstream obj;
				ifstream obj2;
				long double data;
				string filename;
				string filename2;
				filename="final_A_"+to_string(s)+".txt";
				filename2="final_B_"+to_string(s)+".txt";
				obj.open(filename);
				while(obj>>data)
				{
					if(j==N)
					{
						j=0;
						i++;
					}
					A[i][j]=data;
					j++;
				}
				obj.close();

				i=0,j=0;
				obj2.open(filename2);
				while(obj2>>data)
				{
					if(j==M)
					{
						j=0;
						i++;
					}
					B[i][j]=data;
					j++;
				}
				obj2.close();
				forward(obs_seq);
				xander++;
			}

		//DECIDING DIGIT BY FINDING MAX PROBABILITY
		long double the_max=0;
		int the_index=0;
		//2
		for(int d=0;d<8;d++)
		{
			if(d==0)
				cout<<"prob paint     "<<prob_seq[d]<<endl;

			if(d==1)
				cout<<"prob notepad   "<<prob_seq[d]<<endl;

			if(d==2)
				cout<<"prob word      "<<prob_seq[d]<<endl;

			if(d==3)
				cout<<"prob excel     "<<prob_seq[d]<<endl;
 
			if(d==4)
				cout<<"prob photos    "<<prob_seq[d]<<endl;

			if(d==5)
				cout<<"prob skype     "<<prob_seq[d]<<endl;

			//3
			if(d==6)
				cout<<"prob game      "<<prob_seq[d]<<endl;
			
			if(d==7)
				cout<<"prob web       "<<prob_seq[d]<<endl;

			if(d==8)
				cout<<"prob video     "<<prob_seq[d]<<endl;

			if(prob_seq[d]>the_max)
			{
				the_max=prob_seq[d];
				the_index=d;
			}
			
		}

		if(the_max==0)
		{
			//DO NOTHING
		}
		else
		{
			if(the_index==0)
			{
				cout<<"Word is : paint "<<endl;
				STARTUPINFO startInfo={0};
				PROCESS_INFORMATION processInfo={0};
				BOOL bSuccess=CreateProcess(TEXT("C:\\WINDOWS\\system32\\mspaint.exe"),NULL,NULL,NULL,FALSE,NULL,NULL,NULL,&startInfo,&processInfo);
			}
			if(the_index==1)
			{
				cout<<"Word is : notepad "<<endl;
				STARTUPINFO startInfo={0};
				PROCESS_INFORMATION processInfo={0};
				BOOL bSuccess=CreateProcess(TEXT("C:\\WINDOWS\\system32\\notepad.exe"),NULL,NULL,NULL,FALSE,NULL,NULL,NULL,&startInfo,&processInfo);
			}

	
			if(the_index==2)
			{
				cout<<"Word is : word "<<endl;
				STARTUPINFO startInfo={0};
				PROCESS_INFORMATION processInfo={0};
				BOOL bSuccess=CreateProcess(TEXT("C:\\Program Files (x86)\\Microsoft Office\\root\\Office16\\WINWORD.exe"),NULL,NULL,NULL,FALSE,NULL,NULL,NULL,&startInfo,&processInfo);
			}

			if(the_index==3)
			{
				cout<<"Word is : excel "<<endl;
				STARTUPINFO startInfo={0};
				PROCESS_INFORMATION processInfo={0};
				BOOL bSuccess=CreateProcess(TEXT("C:\\Program Files (x86)\\Microsoft Office\\root\\Office16\\EXCEL.exe"),NULL,NULL,NULL,FALSE,NULL,NULL,NULL,&startInfo,&processInfo);
			}
		
			if(the_index==4)
			{
				cout<<"Word is : photos "<<endl;
				STARTUPINFO startInfo={0};
				PROCESS_INFORMATION processInfo={0};
				BOOL bSuccess=CreateProcess(TEXT("C:\\Program Files (x86)\\FastStone Image Viewer\\FSViewer.exe"),NULL,NULL,NULL,FALSE,NULL,NULL,NULL,&startInfo,&processInfo);
			}

			if(the_index==5)
			{
				cout<<"Word is : skype "<<endl;
				STARTUPINFO startInfo={0};
				PROCESS_INFORMATION processInfo={0};
				BOOL bSuccess=CreateProcess(TEXT("C:\\Program Files (x86)\\Microsoft\\Skype for Desktop\\Skype.exe"),NULL,NULL,NULL,FALSE,NULL,NULL,NULL,&startInfo,&processInfo);
			}

			//4
			if(the_index==6)
			{
				cout<<"Word is : game "<<endl;
				STARTUPINFO startInfo={0};
				PROCESS_INFORMATION processInfo={0};
				BOOL bSuccess=CreateProcess(TEXT("F:\\DOWNLOADS\\tictactoe.exe"),NULL,NULL,NULL,FALSE,NULL,NULL,NULL,&startInfo,&processInfo);
			}

			if(the_index==7)
			{
				cout<<"Word is : web "<<endl;
				STARTUPINFO startInfo={0};
				PROCESS_INFORMATION processInfo={0};
				BOOL bSuccess=CreateProcess(TEXT("C:\\Program Files (x86)\\Google\\Chrome\\Application\\chrome.exe"),NULL,NULL,NULL,FALSE,NULL,NULL,NULL,&startInfo,&processInfo);
			}

			if(the_index==8)
			{
				cout<<"Word is : video "<<endl;
				STARTUPINFO startInfo={0};
				PROCESS_INFORMATION processInfo={0};
				BOOL bSuccess=CreateProcess(TEXT("C:\\Program Files\\VideoLAN\\VLC\\vlc.exe"),NULL,NULL,NULL,FALSE,NULL,NULL,NULL,&startInfo,&processInfo);
			}
		}
}

//TO INPUT TEST FILES
//NOT USED DURING LIVE RECORDING
void input_test(int y,int q)
{
	long long int h=y;
	long long int f=q+21;
	long double data1;
	int z=0;
	count_samples=0;
	ifstream obj;

	string filename;
		filename="194101013_"+to_string(h)+"_"+to_string(f)+".txt";
	
	obj.open(filename);
	
	while(obj>>data1)
	{
		initial_samples[z]=data1;
		//cout<<initial_samples[z]<<endl;
		z++;
	}
	count_samples=z;
	//cout<<"number of samples : "<<count_samples<<endl;

}



void live_train()
{
	long double live_samples[100000];

	for(int i=1;i<=5;i++)
	{
		long double data1;
		int z=0;

		if(i==1)
			system("Recording_Module.exe 3 input_file.wav live_input_file_1.txt");
		if(i==2)
			system("Recording_Module.exe 3 input_file.wav live_input_file_2.txt");
		if(i==3)
			system("Recording_Module.exe 3 input_file.wav live_input_file_3.txt");
		if(i==4)
			system("Recording_Module.exe 3 input_file.wav live_input_file_4.txt");
		if(i==5)
			system("Recording_Module.exe 3 input_file.wav live_input_file_5.txt");	
	}


	for(int count=0;count<3;count++)
		{
			//cout<<endl<<"DIGIT "<<j<<" "<<"RUN NO: "<<count+1<<endl;
			initialize_temp_A_B();
			//FOR 20 UTTERANCE OF SAME DIGIT
			for(int i=1;i<=5;i++)
			{
				live_input(i);
				do_dc_shift();
				Normalise();
				select_frames();       
				find_observation_sequence();

				if(count==0)
				{
					 input_model();  //TO INPUT THE MODEL AND OBSERVATION SEQUENCES
				}
				else
				{
					input_model2();
				}
		      
				int turn=0;
				
				while(turn!=20)
				{
					//cout<<"*******************************************************************"<<endl;
					//cout<<"Iteration Number : "<<turn+1<<endl;
					forward(obs_seq);   //TO RUN FORWARD PROCEDURE
					//forward(obs_seq);  //TO RUN FORWARD PROCEDURE
					backward(obs_seq); //TO RUN BACKWARD PROCEDURE
					calc_gamma();    //TO CALCULATE GAMMA
					viterbi(obs_seq);  //TO CALCULATE DELTA AND SI
					reestimation(obs_seq);
					//cout<<endl;
					turn++;
				}

				//KEEP ON ADDING TO TEMP_A AND TEMP_B
				for(int r=0;r<N;r++)
				{
					for(int t=0;t<N;t++)
					{
						temp_A[r][t]=temp_A[r][t]+A[r][t];
						//cout<<temp_A[r][t]<<"  ";
					}
					//cout<<endl;
				}

				for(int r=0;r<N;r++)
				{
					for(int t=0;t<M;t++)
					{
						temp_B[r][t]=temp_B[r][t]+B[r][t];
						//cout<<temp_B[r][t]<<"  ";
					}
					//cout<<endl;
				}
			}

			//DIVIDE TEMP_A AND TEMP_B BY 20
			for(int r=0;r<N;r++)
				{
					for(int t=0;t<N;t++)
					{
						temp_A[r][t]=temp_A[r][t]/5;
						pass_A[r][t]=temp_A[r][t];
						//cout<<temp_A[r][t]<<"  ";
					}
					//cout<<endl;
				}

			for(int r=0;r<N;r++)
				{
					for(int t=0;t<M;t++)
					{
						temp_B[r][t]=temp_B[r][t]/5;
						pass_B[r][t]=temp_B[r][t];
						//cout<<temp_B[r][t]<<"  ";
					}
					//cout<<endl;
				}
		}

		//cout<<endl<<endl;
	
	//STORE FINAL MODELS
	live_move_to_text();  //TO MOVE ALPHA,BETA,GAMMA,DELTA,SI,XI TO TEXT FILE
}



void live_input(int i)
{
	long long int h=i;
	long double data1;
	int z=0;
	count_samples=0;

	ifstream obj;
	string filename;
		filename="live_input_file_"+to_string(h)+".txt";
	
	obj.open(filename);
	
	while(obj>>data1)
	{
		//cout<<data1<<endl;
		initial_samples[z]=data1;
		//cout<<initial_samples[z]<<endl;
		z++;
	}
	count_samples=z;
	//cout<<"number of samples : "<<count_samples<<endl;
}


void live_move_to_text()
{
	//TO MOVE NEW MODEL TO TEXT FILE
	ofstream obj;
	string filename;
		filename="live_final_A_.txt";
	obj.open(filename);

	for(int i=0;i<N;i++)
	{
		for(int j=0;j<N;j++)
		{
			obj<<temp_A[i][j]<<"    ";
		}
		obj<<endl;
	}

	ofstream obj2;
	string filename2;
		filename2="live_final_B_.txt";
	obj2.open(filename2);

	for(int i=0;i<N;i++)
	{
		for(int j=0;j<M;j++)
		{
			obj2<<temp_B[i][j]<<"    ";
		}
		obj2<<endl;
	}
}