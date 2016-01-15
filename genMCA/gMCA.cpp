#include<iostream>
#include<unordered_map>
#include<vector>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<algorithm>
#include<windows.h>
using namespace std;
bool findMCA;
int p;
int t,N,k,paraCombNum,curCombNum,totalCombNum,minCombRow,minRowIndex;
double T;
int existCombNum;
vector<vector<int>> MCA;
vector<vector<bool>> combinations;
vector<int*> combIndex;
vector<int> combNumInRow;
vector<vector<bool>> combWhetherUnique;
vector<int> v;
vector<int> runtime;

class key
{
public:
	key(int* p1,int p2)
	{
		value=new int[t];
		for(int i=0;i<t;i++)
		{
			value[i]=p1[i];
		}
		paraNum=p2;
	}

	bool operator == (const key& v1)
	{
		if(v1.paraNum!=paraNum)
			return false;
		else
		{
			for(int i=0;i<t;i++)
			{
				if(v1.value[i]!=value[i])
					return false;
			}
			return true;
		}
	}

	int* value;
	int paraNum;
};
class hash_class
{
public:
	size_t operator() (const key& k) const
    {
		int hashcode=k.value[0];
		for(int i=1;i<t;i++)
			hashcode=hashcode*v[i]+k.value[i];
		return (size_t)hashcode;
    }
};
class equal_class
{
public:
	bool operator() (const key& v1,const key& v2) const
	{
		if(v1.paraNum!=v2.paraNum)
			return false;
		else
		{
			for(int i=0;i<t;i++)
			{
				if(v1.value[i]!=v2.value[i])
					return false;
			}
			return true;
		}
	}
};

vector<vector<key>> combInRow;
unordered_map <key,int,hash_class,equal_class> G;

int calculateHD(int row,vector<int> curCandidateRow)
{
	int hmDistance=0;
	for(int j=0;j<k;j++)
	{
		int columnHMD=0;
		for(int i=0;i<row;i++)
		{
			if(MCA[i][j]!=curCandidateRow[j])
				columnHMD++;
		}
		hmDistance=hmDistance+columnHMD;
	}
	return hmDistance;
}
void initMCA()
{
	srand(time(NULL));
	vector<int> firstRow;
	for(int j=0;j<k;j++)
		firstRow.insert(firstRow.begin()+j,rand()%v[j]);
	MCA.insert(MCA.begin(),firstRow);
	for(int row=1;row<N;row++)
	{
		vector<int> firstTry;
		for(int j=0;j<k;j++)
			firstTry.insert(firstTry.begin()+j,rand()%v[j]);
		MCA.insert(MCA.begin()+row,firstTry);
		int maxHDRow=calculateHD(row,firstTry);
		for(int i=1;i<p;i++)
		{
			vector<int> temp;
			for(int j=0;j<k;j++)
				temp.insert(temp.begin()+j,rand()%v[j]);
			int curRowHD=calculateHD(row,temp);
			if(curRowHD>maxHDRow)
			{
				MCA[row]=temp;
				maxHDRow=curRowHD;
			}
		}
		//printf("maxHammingDistance is %d\n",maxHDRow);
	}
}
void getCombIndex(int c,int* p)
{
	int* index=new int[t];
	for(int i=0,j=0;i<k;i++)
	{
		if(p[i]==1)
			index[j++]=i;
	}
	combIndex.insert(combIndex.begin()+c,index);
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
		getCombIndex(paraCombNum,perIndex);
		paraCombNum++;
	}while(next_permutation(perIndex,perIndex+k));
}
void addCombToG(vector<int> newRow)
{
	for(int c=0;c<paraCombNum;c++)
	{
		int* value=new int[t];
		for(int i=0;i<t;i++)
		{
			value[i]=newRow[combIndex[c][i]];
		}
		key comb=key(value,c);
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
		vector<key> temp2;
		int* value=new int[t];
		for(int i=0;i<t;i++)
		{
			value[i]=MCA[0][i];
		}
		for(int j=0;j<paraCombNum;j++)
		{
			temp2.insert(temp2.begin()+j,key(value,-1));
		}
		combInRow.insert(combInRow.begin()+i,temp2);
		combNumInRow.insert(combNumInRow.begin()+i,0);
	}
}
void decideCIR()
{
	minCombRow=paraCombNum;
	minRowIndex=-1;
	for(int row=0;row<N;row++)
	{
		int tempMinCombRow=0; 
		for(int c=0;c<paraCombNum;c++)
		{
			int* value=new int[t];
			for(int i=0;i<t;i++)
			{
				value[i]=MCA[row][combIndex[c][i]];
			}
			key comb=key(value,c);
			if(G.find(comb)!=G.end())
			{
				combInRow[row][c]=comb;
				if(G[comb]==1)
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
void updateCIR(int type,vector<int> row,int curRowNum,int rowNum,bool ifNeed)// type=1 erase old row,=2 replace old row,=3 add new row;curRowNum is the parameter row's num,rowNum is new row's num
{
	if(type==1)
	{
		combNumInRow[curRowNum]=0;
	}
	else if(type==3)
	{
		combNumInRow[rowNum]=0;
	}
	for(int c=0;c<paraCombNum;c++)
	{
		int* value=new int[t];
		for(int i=0;i<t;i++)
		{
			value[i]=row[combIndex[c][i]];
		}
		key comb=key(value,c);
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
						if(combInRow[i][c]==comb&&i!=curRowNum)
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
			combInRow[curRowNum][c]=comb;
			if(G.find(comb)==G.end())
			{
				G.insert(make_pair(comb,1));
				curCombNum++;
				combWhetherUnique[curRowNum][c]=true;
				combNumInRow[curRowNum]++;
			}
			else
			{
				G[comb]++;
				combWhetherUnique[curRowNum][c]=false;
				if(G[comb]==2)
				{
					for(int i=0;i<N;i++)
					{
						if(combInRow[i][c]==comb&&i!=curRowNum)
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
						if(combInRow[i][c]==comb)
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
	if(ifNeed==false)
		return;
	minCombRow=paraCombNum;
	int r=N;
	if(type==3)
		r=r+1;
	for(int i=0;i<r;i++)
	{
		int temp=combNumInRow[i];
		if(temp<minCombRow)
		{
			minCombRow=temp;
			minRowIndex=i;
		}
	}
}
int fitness(int type,int oldRowIndex,vector<int> individual)//type=1 judge if we get a MCA,=2 the fitness with the old row,=3 the fitness after replace the old row with the new row
{
	if(type==1||type==2)
	{
		return totalCombNum-curCombNum;
	}
	else
	{
		int incresedComb=0; 
		for(int c=0;c<paraCombNum;c++)
		{
			int* value=new int[t];
			for(int i=0;i<t;i++)
			{
				value[i]=individual[combIndex[c][i]];
			}
			key comb=key(value,c);
			if(G.find(comb)==G.end()||(combInRow[oldRowIndex][c]==comb&&combWhetherUnique[oldRowIndex][c]==true))
			{
				incresedComb++;
			}
		}
		return totalCombNum-(curCombNum-combNumInRow[oldRowIndex])-incresedComb;
	}
}
void addOneRow(int rowNum)
{
	vector<int> firstTry;
	for(int j=0;j<k;j++)
		firstTry.insert(firstTry.begin()+j,rand()%v[j]);
	MCA.insert(MCA.begin()+rowNum,firstTry);
	int maxHDRow=calculateHD(rowNum,firstTry);
	for(int i=1;i<p;i++)
	{
		vector<int> temp;
		for(int j=0;j<k;j++)
			temp.insert(temp.begin()+j,rand()%v[j]);
		int curRowHD=calculateHD(rowNum,temp);
		if(curRowHD>maxHDRow)
		{
			MCA[rowNum]=temp;
			maxHDRow=curRowHD;
		}
	}
	vector<bool> temp1;
	for(int j=0;j<paraCombNum;j++)
	{
		temp1.insert(temp1.begin()+j,false);
	}
	combWhetherUnique.insert(combWhetherUnique.begin()+rowNum,temp1);
	vector<key> temp2;
	int* value=new int[t];
	for(int i=0;i<t;i++)
	{
		value[i]=MCA[0][i];
	}
	for(int j=0;j<paraCombNum;j++)
	{
		temp2.insert(temp2.begin()+j,key(value,-1));
	}
	combInRow.insert(combInRow.begin()+rowNum,temp2);
	combNumInRow.insert(combNumInRow.begin()+rowNum,0);
}
void deleteZRow()
{
	for(int i=0;i<N;i++)
	{
		if(combNumInRow[i]==0)
		{
			updateCIR(1,MCA[i],i,0,false);
			std::vector<vector<int>>::iterator iter1=MCA.begin()+i;
			MCA.erase(iter1);
			std::vector<vector<key>>::iterator iter2=combInRow.begin()+i;
			combInRow.erase(iter2);
		}
	}
}
void printMCA()
{
	deleteZRow();
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<k;j++)
		{
			printf("%d ",MCA[i][j]);
		}
		printf("\n");
	}
	printf("totalCombNum is %d,N is %d.\n",totalCombNum,N);
}
double generateP(int c1,int c2,int T)
{
	return exp(((double)(c1-c2))/((double)T));
}
void dfs(vector<int> newRow,int* curComb,int* vNum,int dep,int c,bool& flag)
{
	if(flag==true)
		return;
	if(dep<t)
	{
		for(int i=0;i<vNum[dep];i++)
		{
			curComb[dep]=i;
			dfs(newRow,curComb,vNum,dep+1,c,flag);
			if(flag==true)
				return;
		}
	}
	else
	{
		key comb=key(curComb,c);
		if(G.find(comb)==G.end())
		{
			for(int i=0;i<t;i++)
			{
				newRow[combIndex[c][i]]=curComb[i];
			}
			int rowIndex,curRes=0,minRes=100000000;
			int c1=fitness(2,0,MCA[0]);
			int c2;
			for(int row=0;row<N;row++)
			{
				for(int i=0;i<k;i++)
				{
					if(combinations[c][i]==false)
					{
						newRow[i]=MCA[row][i];
					}
				}
				c2=fitness(3,row,newRow);
				int res=c2-c1;
				if(res<minRes)
				{
					minRes=res;
					rowIndex=row;
				}
			}
			double tempp=generateP(c1,minRes+c1,T);
			if(minRes<curRes||(tempp>(double)(rand()%10000)/10000.0))
			{
				for(int i=0;i<k;i++)
				{
					if(combinations[c][i]==false)
					{
						newRow[i]=MCA[rowIndex][i];
					}
				}
				vector<int> temp=MCA[rowIndex];
				MCA[rowIndex]=newRow;
				updateCIR(1,temp,rowIndex,0,false);
				updateCIR(2,newRow,rowIndex,0,false);
				flag=true;
			}
		}
		else
		{
			existCombNum++;
		}
	}
	return;
}
void tryAddOneTuple(int c)
{
	int* curComb=new int[t];
	int* vNum=new int[t];
	int totalComb=1;
	for(int i=0;i<t;i++)
	{
		vNum[i]=v[combIndex[c][i]];
		totalComb=totalComb*vNum[i];
	}
	bool flag=false;
	existCombNum=0;
	vector<int> newRow;
	for(int i=0;i<k;i++)
	{
		newRow.insert(newRow.begin()+i,0);
	}
	dfs(newRow,curComb,vNum,0,c,flag);
	while(existCombNum==totalComb&&fitness(1,0,MCA[0])!=0)
	{
		c=(c+1)%paraCombNum;
		int totalComb=1;
		for(int i=0;i<t;i++)
		{
			vNum[i]=v[combIndex[c][i]];
			totalComb=totalComb*vNum[i];
		}
		bool flag=false;
		existCombNum=0;
		dfs(newRow,curComb,vNum,0,c,flag);
	}
	return;
}
void tryChangeMij()
{
	int row=rand()%N;
	int column=rand()%k;
	int curMij=MCA[row][column];
	vector<int> newRow=MCA[row];
	int curRes=0,minRes=100000000;
	int maxColumnv=curMij;
	int c1=fitness(2,row,MCA[0]);
	int c2;
	for(int i=0;i<v[column];i++)
	{
		if(i!=curMij)
		{
			newRow[column]=i;
			c2=fitness(3,row,newRow);
			int res=c2-c1;
			if(res<minRes)
			{
				minRes=res;
				maxColumnv=i;
			}
		}
	}
	int tempp=generateP(c1,minRes+c1,T);
	if(minRes<curRes||(tempp>((double)(rand()%10000))/10000.0))
	{
		newRow[column]=maxColumnv;
		vector<int> temp=MCA[row];
		MCA[row]=newRow;
		updateCIR(1,temp,row,0,false);
		updateCIR(2,newRow,row,0,false);
	}
}
void verifyMCA()
{
	G.clear();
	for(int i=0;i<N;i++)
	{
		addCombToG(MCA[i]);
	}
	if(G.size()==totalCombNum)
		printf("The answer is correct!\n");
	else
		printf("The answer is wrong!\n");
}
void getRuntime()
{
	int allPossibleTestCaseNum=1;
	for(int i=0;i<k;i++)
		allPossibleTestCaseNum=allPossibleTestCaseNum*v[i];
	for(int j=0;j<allPossibleTestCaseNum;j++)
		runtime.insert(runtime.begin()+j,rand());
}
int getRowIndex(int rowNum)
{
	int index=MCA[rowNum][0];
	for(int i=1;i<k;i++)
	{
		index=index*v[i]+MCA[rowNum][i];
	}
	return index;
}
int getMaxRowIndex()
{
	int index=0;
	for(int i=1;i<N;i++)
	{
		if(runtime[getRowIndex(i)]>runtime[getRowIndex(index)])
			index=i;
	}
	return index;
}
void replaceOneRow()
{
	//find all combs that are unique, the results can be used to construct one row(or two rows, each has half of the unique combs) that has all the unique combs in the old row while its runtime is less.

}
int main()
{
	//initialize parameters and MCA
	totalCombNum=0;
	curCombNum=0;
	paraCombNum=0;
	printf("Input parameter t k p v1 v2 ... vk:\n");
	scanf("%d%d%d",&t,&k,&p);
	if(k==1)
	{
		printf("error!\n");
		return 0;
	}
	int* vt=new int[k];
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
		vt[i]=temp;
	}
	for(int i=1;i<k;i++)
	{
		for(int j=0;j<i;j++)
		{
			if(vt[j]<vt[i])
			{
				int temp=vt[j];
				vt[j]=vt[i];
				vt[i]=temp;
			}
		}
	}
	N=1;
	for(int i=0;i<t;i++)
	{
		N=N*vt[i];
	}
	initMCA();
	genCombs();
	getRuntime();
	for(int i=0;i<N;i++)
	{
		addCombToG(MCA[i]);
	}
	if(fitness(1,0,MCA[0])==0)
	{
		printMCA();
		verifyMCA();
		return 0;
	}
	findMCA=false;
	initCIR();
	decideCIR();
	//Simulated Annealing algorithm
	double T0=4.0,Tf=pow(10,-10),coolingFactor=0.99;
	int V=N,frozenFactor=11;
	int L;
	int tDecreWOCombNumChange;
	T=T0;
	while(findMCA!=true)
	{
		tDecreWOCombNumChange=0;
		int bestSFCurNum=curCombNum;
		while(findMCA!=true&&tDecreWOCombNumChange<frozenFactor&&T>Tf)
		{
			//L=N*k*V*V;
			L=N*k*V;
			int storeCurNum=bestSFCurNum;
			//int count=0,l=N*k*V;
			for(int i=0;i<L;i++)
			{
				double r=((double)(rand()%10000))/10000.0;
				if(r<0.7)
				{
					tryChangeMij();
				}
				else
				{
					int c=rand()%paraCombNum;
					tryAddOneTuple(c);
				}
				if(fitness(1,0,MCA[0])==0)
				{
					printMCA();
					verifyMCA();
					return 0;
				}
				else
				{
					if(curCombNum>bestSFCurNum)
					{
						bestSFCurNum=curCombNum;
					}
					/*else
					{
						count++;
					}
					if(count>l)
						break;*/
				}
			}
			T=T*coolingFactor;
			if(bestSFCurNum==storeCurNum)
				tDecreWOCombNumChange++;
			else
				tDecreWOCombNumChange=0;
		}
		addOneRow(N);
		updateCIR(3,MCA[N],0,N,true);
		N=N+1;
		if(fitness(1,0,MCA[0])==0)
		{
			printMCA();
			verifyMCA();
			return 0;
		}
	}
	return 0;
}