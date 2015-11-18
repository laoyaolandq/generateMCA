#include<iostream>
#include<hash_map>
#include<vector>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<algorithm>
#include<windows.h>
#include<StringAsKey.h>
#include<algorithm>
using namespace std;
bool findMCA;
int t,N,k,maxV,secMaxV,paraCombNum,curCombNum,totalCombNum,minCombRow,minRowIndex;
int NP,gNum;
float F,CR;
vector<vector<int>> MCA;
vector<vector<int>> group;
vector<vector<bool>> combinations;
vector<vector<string>> combInRow;
vector<int> combNumInRow;
vector<vector<bool>> combWhetherUnique;
vector<int> v;
hash_map <string,int> G;
void initMCA()
{
	srand(time(NULL)+rand());
	for(int i=0;i<N;i++)
	{
		vector<int> temp;
		for(int j=0;j<k;j++)
			temp.insert(temp.begin()+j,rand()%v[j]);
		MCA.insert(MCA.begin()+i,temp);
	}
}
void initGroup()
{
	srand(time(NULL)+rand());
	for(int i=0;i<NP;i++)
	{
		vector<int> temp;
		for(int j=0;j<k;j++)
			temp.insert(temp.begin()+j,rand()%v[j]);
		group.insert(group.begin()+i,temp);
	}
}
void genCombs()
{
	int* perIndex=new int[k];
	for(int i=0;i<k;i++)
	{
		if(i<t)
			perIndex[i]=1;
		else
			perIndex[i]=0;
	}
	reverse(perIndex,perIndex+k);
	do
	{
		vector<bool> temp;
		int tempCombNum=1;
		for(int i=0;i<k;i++)
		{
			if(perIndex[i]==1)
			{
				temp.insert(temp.begin()+i,true);
				tempCombNum*=v[i];
			}
			else
				temp.insert(temp.begin()+i,false);
		}
		totalCombNum+=tempCombNum;
		combinations.insert(combinations.begin()+paraCombNum,temp);
		paraCombNum++;
	}while(next_permutation(perIndex,perIndex+k));
}
void addCombToG(vector<int> newRow)
{
	for(int c=0;c<paraCombNum;c++)
	{
		string comb;
		int curIndex=0;
		for(int i=k-1;i>=0;i--)
		{
			if(combinations[c][i]==false)
			{
				int tempV,zeroNum=0;
				tempV=v[i]-1;
				do
				{
					zeroNum++;
					tempV=tempV/10;
				}while(tempV!=0);
				comb.insert(curIndex,zeroNum,'x');
				curIndex+=zeroNum;
			}
			else
			{
				int tempV,tempMij;
				tempV=v[i]-1;
				tempMij=newRow[i];
				do
				{
					tempV=tempV/10;
					int temp=tempMij%10;
					comb.insert(curIndex++,1,(char)(temp+48));
					tempMij=tempMij/10;
				}while(tempV!=0);
			}
		}
		if(G.find(comb)==G.end())
		{
			G.insert(make_pair(comb,1));
			curCombNum++;
		}
		else
			G[comb]++;
	}
}
void initCIR()
{
	for(int i=0;i<N;i++)
	{
		vector<bool> temp1;
		for(int j=0;j<paraCombNum;j++)
		{
			temp1.insert(temp1.begin()+j,false);
		}
		combWhetherUnique.insert(combWhetherUnique.begin()+i,temp1);
		vector<string> temp2;
		for(int j=0;j<paraCombNum;j++)
		{
			temp2.insert(temp2.begin()+j,"");
		}
		combInRow.insert(combInRow.begin()+i,temp2);
		combNumInRow.insert(combNumInRow.begin()+i,0);
	}
}
void decideCIR()
{
	for(int row=0;row<N;row++)
	{
		int tempMinCombRow=0; 
		for(int c=0;c<paraCombNum;c++)
		{
			string comb;
			int curIndex=0;
			for(int i=k-1;i>=0;i--)
			{
				if(combinations[c][i]==false)
				{
					int tempV,zeroNum=0;
					tempV=v[i]-1;
					do
					{
						zeroNum++;
						tempV=tempV/10;
					}while(tempV!=0);
					comb.insert(curIndex,zeroNum,'x');
					curIndex+=zeroNum;
				}
				else
				{
					int tempV,tempMij;
					tempV=v[i]-1;
					tempMij=MCA[row][i];
					do
					{
						tempV=tempV/10;
						int temp=tempMij%10;
						comb.insert(curIndex++,1,(char)(temp+48));
						tempMij=tempMij/10;
					}while(tempV!=0);
				}
			}
			hash_map<string, int>::iterator iter=G.find(comb);
			if(iter!=G.end())
			{
				combInRow[row][c]=comb;
				int valueNum=iter->second;
				if(valueNum==1)
				{
					combWhetherUnique[row][c]=true;
					tempMinCombRow++;
				}
				else
					combWhetherUnique[row][c]=false;
			}
			else
				printf("error\n");
		}
		combNumInRow[row]=tempMinCombRow;
		if(tempMinCombRow<minCombRow)
		{
			minCombRow=tempMinCombRow;
			minRowIndex=row;
		}
	}
}
int countT(vector<bool> w)
{
	int count=0;
	for(int i=0;i<paraCombNum;i++)
		if(w[i]==true)
			count++;
	return count;
}
void updateCIR(int type,vector<int> row,int rowNum)
{
	combNumInRow[minRowIndex]=0;
	for(int c=0;c<paraCombNum;c++)
	{
		string comb;
		int curIndex=0;
		for(int i=k-1;i>=0;i--)
		{
			if(combinations[c][i]==false)
			{
				int tempV,zeroNum=0;
				tempV=v[i]-1;
				do
				{
					zeroNum++;
					tempV=tempV/10;
				}while(tempV!=0);
				comb.insert(curIndex,zeroNum,'x');
				curIndex+=zeroNum;
			}
			else
			{
				int tempV,tempMij;
				tempV=v[i]-1;
				tempMij=row[i];
				do
				{
					tempV=tempV/10;
					int temp=tempMij%10;
					comb.insert(curIndex++,1,(char)(temp+48));
					tempMij=tempMij/10;
				}while(tempV!=0);
			}
		}
		if(type==1)
		{
			if(G[comb]==1)
			{
				G.erase(comb);
				curCombNum--;
			}
			else
			{
				G[comb]--;
				if(G[comb]==1)
				{
					for(int i=0;i<N;i++)
					{
						if(combInRow[i][c]==comb&&i!=minRowIndex)
						{
							combWhetherUnique[i][c]=true;
							combNumInRow[i]++;
							break;
						}
					}
				}
			}
		}
		else if(type==2)
		{
			combInRow[minRowIndex][c]=comb;
			if(G.find(comb)==G.end())
			{
				G.insert(make_pair(comb,1));
				curCombNum++;
				combWhetherUnique[minRowIndex][c]=true;
				combNumInRow[minRowIndex]++;
			}
			else
			{
				G[comb]++;
				combWhetherUnique[minRowIndex][c]=false;
				if(G[comb]==2)
				{
					for(int i=0;i<N;i++)
					{
						if(combInRow[i][c]==comb&&i!=minRowIndex)
						{
							combWhetherUnique[i][c]=false;
							combNumInRow[i]--;
							break;
						}
					}
				}
			}
		}
		else
		{
			combInRow[rowNum][c]=comb;
			if(G.find(comb)==G.end())
			{
				G.insert(make_pair(comb,1));
				curCombNum++;
				combWhetherUnique[rowNum][c]=true;
				combNumInRow[rowNum]++;
			}
			else
			{
				G[comb]++;
				combWhetherUnique[rowNum][c]=false;
				if(G[comb]==2)
				{
					for(int i=0;i<N;i++)
					{
						if(combInRow[i][c]==comb&&i!=rowNum)
						{
							combWhetherUnique[i][c]=false;
							combNumInRow[i]--;
							break;
						}
					}
				}
			}
		}
	}
	minCombRow=paraCombNum;
	for(int i=0;i<N;i++)
	{
		int temp=combNumInRow[i];
		if(temp<minCombRow)
		{
			minCombRow=temp;
			minRowIndex=i;
		}
		/*int temp=countT(combWhetherUnique[i]);
		combNumInRow[i]=temp;
		if(temp<minCombRow)
		{
			minCombRow=temp;
			minRowIndex=i;
		}*/
	}
}
int fitness(int type,vector<int> individual)
{
	int missedCom=totalCombNum-curCombNum;
	if(type==1||type==2)
	{
		return missedCom;
	}
	else
	{
		int incresedComb=0; 
		for(int c=0;c<paraCombNum;c++)
		{
			string comb;
			int curIndex=0;
			for(int i=k-1;i>=0;i--)
			{
				if(combinations[c][i]==false)
				{
					int tempV,zeroNum=0;
					tempV=v[i]-1;
					do
					{
						zeroNum++;
						tempV=tempV/10;
					}while(tempV!=0);
					comb.insert(curIndex,zeroNum,'x');
					curIndex+=zeroNum;
				}
				else
				{
					int tempV,tempMij;
					tempV=v[i]-1;
					tempMij=individual[i];
					do
					{
						tempV=tempV/10;
						int temp=tempMij%10;
						comb.insert(curIndex++,1,(char)(temp+48));
						tempMij=tempMij/10;
					}while(tempV!=0);
				}
			}
			if(G.find(comb)==G.end()||(combInRow[minRowIndex][c]==comb&&combWhetherUnique[minRowIndex][c]==true))
			{
				incresedComb++;
			}
		}
		return totalCombNum-(curCombNum-minCombRow)-incresedComb;
	}
}
void addOneRow(int rowNum)
{
	vector<int> temp;
	for(int j=0;j<k;j++)
		temp.insert(temp.begin()+j,rand()%v[j]);
	MCA.insert(MCA.begin()+rowNum,temp);
	vector<bool> temp1;
	for(int j=0;j<paraCombNum;j++)
	{
		temp1.insert(temp1.begin()+j,false);
	}
	combWhetherUnique.insert(combWhetherUnique.begin()+rowNum,temp1);
	vector<string> temp2;
	for(int j=0;j<paraCombNum;j++)
	{
		temp2.insert(temp2.begin()+j,"");
	}
	combInRow.insert(combInRow.begin()+rowNum,temp2);
	combNumInRow.insert(combNumInRow.begin()+rowNum,0);
}
void printMCA()
{
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<k;j++)
		{
			printf("%d ",MCA[i][j]);
		}
		printf("\n");
	}
}
void addOneTuple(int c)
{
	hash_map<string, int> specTuple;
	int existComNum=0,maxNum=-1;
	string maxNumIndex;
	for(int row=0;row<N;row++)
	{
		string comb;
		int curIndex=0;
		for(int i=k-1;i>=0;i--)
		{
			if(combinations[c][i]==false)
			{
				int tempV,zeroNum=0;
				tempV=v[i]-1;
				do
				{
					zeroNum++;
					tempV=tempV/10;
				}while(tempV!=0);
				comb.insert(curIndex,zeroNum,'x');
				curIndex+=zeroNum;
			}
			else
			{
				int tempV,tempMij;
				tempV=v[i]-1;
				tempMij=MCA[row][i];
				do
				{
					tempV=tempV/10;
					int temp=tempMij%10;
					comb.insert(curIndex++,1,(char)(temp+48));
					tempMij=tempMij/10;
				}while(tempV!=0);
			}
		}
		hash_map<string, int>::iterator iter=specTuple.find(comb);
		if(iter==specTuple.end())
		{
			specTuple.insert(make_pair(comb,1));
			existComNum++;
		}
		else
		{
			specTuple[comb]++;
			if(specTuple[comb]>maxNum)
			{
				maxNum=specTuple[comb];
				maxNumIndex=comb;
			}
		}
	}
	int total=1;
	for(int i=k-1;i>=0;i--)
	{
		if(combinations[c][i]==true)
		{
			total=total*v[i];
		}
	}
	if(existComNum==total)
		return;

}
int main()
{
	F=0.2;
	CR=0.2;
	totalCombNum=0;
	curCombNum=0;
	paraCombNum=0;
	maxV=-1;
	secMaxV=-1;
	printf("Input parameter t k NP GNUM v1 v2 ... vk:\n");
	scanf("%d%d%d%d",&t,&k,&NP,&gNum);
	if(k==1)
	{
		printf("error!\n");
		return 0;
	}
	/*int a[4]={0,0,1,1};
	do
	{
		printf("%d%d%d%d",a[0],a[1],a[2],a[3]);
	}while(next_permutation(a,a+4));*/
	for(int i=0;i<k;i++)
	{
		int temp;
		scanf("%d",&temp);
		if(temp==0)
		{
			printf("error!\n");
			return 0;
		}
		v.insert(v.begin()+i,temp);
		if(v[i]>maxV)
		{
			secMaxV=maxV;
			maxV=v[i];
			continue;
		}
		if(v[i]==maxV&&v[i]>secMaxV)
		{
			secMaxV=v[i];
			continue;
		}
		if(v[i]<maxV&&v[i]>secMaxV)
		{
			secMaxV=v[i];
			continue;
		}
	}
	//initialize parameters and MCA
	N=maxV*secMaxV;
	LARGE_INTEGER  large_interger;  
	double dff;  
	__int64  c1, c2;  
	QueryPerformanceFrequency(&large_interger);  
	dff = large_interger.QuadPart;  
	QueryPerformanceCounter(&large_interger);  
	c1 = large_interger.QuadPart;
	initMCA();
	genCombs();
	for(int i=0;i<N;i++)
	{
		addCombToG(MCA[i]);
	}
	if(fitness(1,MCA[0])==0)
	{
		printMCA();
		return 0;
	}
	findMCA=false;
	minCombRow=paraCombNum;
	minRowIndex=-1;
	initCIR();
	decideCIR();
	while(findMCA!=true)
	{
		int i=0;
		initGroup();//initialize group
		while(i++<gNum)
		{
			for(int crow=0;crow<NP;crow++)
			{
				//mutation
				vector<int> candidateIn,fCandiIn;
				srand(time(NULL)+rand());
				int r1,r2,r3;
				do
				{
					r1=rand()%NP;
				}while(r1==crow);
				do
				{
					r2=rand()%NP;
				}while(r2==crow||r2==r1);
				do
				{
					r3=rand()%NP;
				}while(r3==crow||r3==r1||r3==r2);
				//float t=exp(1-gNum/(gNum+1-(i+1)));
				//float f=F*pow(2,t);
				for(int i=0;i<k;i++)
				{
					int temp=abs((int)(group[r1][i]+F*(group[r2][i]-group[r3][i])))%v[i];
					candidateIn.insert(candidateIn.begin()+i,temp);
				}
				//crossover
				srand(time(NULL)+rand());
				int jrand=rand()%k;
				for(int i=0;i<k;i++)
				{
					if(rand()%10<(CR*10)||i==jrand)
						fCandiIn.insert(fCandiIn.begin()+i,candidateIn[i]);
					else
						fCandiIn.insert(fCandiIn.begin()+i,group[crow][i]);
				}
				//selection
				//vector<int> fCandiIn=group[crow];
				int res=fitness(3,fCandiIn)-fitness(2,MCA[minRowIndex]);
				if(res<0||(res==0&&MCA[minRowIndex]!=fCandiIn))
				{
					vector<int> temp=MCA[minRowIndex];
					MCA[minRowIndex]=fCandiIn;
					updateCIR(1,temp,0);
					updateCIR(2,fCandiIn,0);
				}
				if(fitness(1,MCA[0])==0)
				{
					QueryPerformanceCounter(&large_interger);  
					c2 = large_interger.QuadPart;
					printMCA();
					printf("N=%d ”√ ±%f ms\n",N,(c2 - c1) * 1000 / dff);
					return 0;
				}
			}
		}
		addOneRow(N);
		updateCIR(3,MCA[N],N);
		N=N+1;
	}
	return 0;
}