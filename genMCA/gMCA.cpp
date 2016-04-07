#include<iostream>
#include<unordered_map>
#include<vector>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<algorithm>
#include<windows.h>
#include<fstream>
#include<sstream>
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
vector<int*> cIndex;
vector<int> lessRTR;
vector<vector<int>> group;
vector<int> groupPro;
vector<vector<int>> candidateRows;
int candidateRsRT;
vector<int> cMinRT;
int cMinRTV;
int curRT;
int sumRT;
bool changeFlag;
int NP,gNum;
float F,CR;

class key
{
public:
	key(int* p1,int p2)
	{
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

	int value[6];
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
void addCIndex(int c,int* perIndex)
{
	int* temp=new int[t];
	int cNum=0;
	for(int i=0;i<k;i++)
	{
		if(perIndex[i]==1)
			temp[cNum++]=i;
	}
	cIndex.insert(cIndex.begin()+c,temp);
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
		addCIndex(paraCombNum,perIndex);
		paraCombNum++;
	}while(next_permutation(perIndex,perIndex+k));
	delete[] perIndex;
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
		delete[] value;
		if(G.find(comb)==G.end())
		{
			G.insert(make_pair(comb,1));
			curCombNum++;
		}
		else
			G[comb]++;
		//delete[] comb.value;
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
		delete[] value;
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
			delete[] value;
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
		delete[] value;
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
			delete[] value;
			if(G.find(comb)==G.end()||(combInRow[oldRowIndex][c]==comb&&combWhetherUnique[oldRowIndex][c]==true))
			{
				incresedComb++;
			}
		}
		return totalCombNum-(curCombNum-combNumInRow[oldRowIndex])-incresedComb;
	}
}
void addOneRow(int rowNum,int type,vector<int> newRow)
{
	if(type==1)
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
	}
	else
	{
		MCA.insert(MCA.begin()+rowNum,newRow);
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
	delete[] value;
}
bool deleteZRow()
{
	int flag=0;
	for(int i=0;i<N;i++)
	{
		if(combNumInRow[i]==0)
		{
			flag=1;
			updateCIR(1,MCA[i],i,0,false);
			std::vector<vector<int>>::iterator iter1=MCA.begin()+i;
			MCA.erase(iter1);
			std::vector<vector<key>>::iterator iter2=combInRow.begin()+i;
			combInRow.erase(iter2);
			std::vector<vector<bool>>::iterator iter3=combWhetherUnique.begin()+i;
			combWhetherUnique.erase(iter3);
			std::vector<int>::iterator iter4=combNumInRow.begin()+i;
			combNumInRow.erase(iter4);
			cout<<i<<" is zero."<<endl;
			N=N-1;
		}
	}
	if(flag==1)
		return true;
	else
		return false;
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
	delete[] curComb;
	delete[] vNum;
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
	curCombNum=0;
	for(int i=0;i<N;i++)
	{
		addCombToG(MCA[i]);
	}
	if(G.size()==totalCombNum&&curCombNum==totalCombNum)
		printf("The answer is correct!\n");
	else
		printf("The answer is wrong!\n");
}
void getAllCaseRuntime(int type)
{
	/*long long allPossibleTestCaseNum=1;
	for(int i=0;i<k;i++)
		allPossibleTestCaseNum=allPossibleTestCaseNum*v[i];
	if(type==1)
	{
		for(long long j=0;j<allPossibleTestCaseNum;j++)
			runtime.insert(runtime.begin()+j,rand());
	}
	else
	{
		for(long long j=0;j<allPossibleTestCaseNum;j++)
			runtime.insert(runtime.begin()+j,-1);
	}*/
}
int getRuntime(vector<int> row)
{
	//the function which is used to calculate the rows' runtime is defined here
	/*int RTime=0;
	for(int i=0;i<k;i++)
	{
		RTime=RTime+row[i]*row[i]*row[i];
	}
	return RTime;*/
	int RTime=0;
	/*for(int i=0;i<k-1;i++)
	{
		RTime=RTime+row[i]*row[i+1];
	}*/
	double temp=349.308*row[0]*row[1]+388.4343*row[1]*row[2]+494.8463*row[2]*row[3]+467.3195*row[3]*row[4]+569.1318*row[4]*row[5]+541.4659*row[5]*row[6]+1077.6573*row[6]*row[7]+414.3779*row[7]*row[8]-29.5446*row[8]*row[9]-457.2587*row[9]*row[10];
	RTime=(int)temp;
	return RTime;
	/*long long product=row[0];
	for(int i=1;i<k;i++)
	{
		product=product*v[i]+row[i];
	}
	if(runtime[product]==-1)
	{
		int RTime=0;
		for(int i=0;i<k-1;i++)
		{
			RTime=RTime+row[i]*row[i]*row[i]*row[i+1]*row[i+1]*row[i+1];
		}
		runtime[product]=RTime;
		return RTime;
	}
	else
	{
		return runtime[product];
	}*/
}
void dfs2(vector<int> newRow,int* index,int dep)
{
	if(dep<k)
	{
		if(index[dep]==1)
			dfs2(newRow,index,dep+1);
		else
		{
			for(int i=0;i<v[dep];i++)
			{
				newRow[dep]=i;
				dfs2(newRow,index,dep+1);
			}
		}
	}
	else
	{
		
		int temp=getRuntime(newRow);
		if(temp<curRT)
		{
			changeFlag=true;
			curRT=temp;
			lessRTR=newRow;
		}
	}
}
void replaceOneRow(int rowNum)
{
	//find all combs that are unique, the results can be used to construct one row(or two rows, each has half of the unique combs) that has all the unique combs in the old row while its runtime is less.
	int* index=new int[k];
	for(int i=0;i<k;i++)
	{
		index[i]=0;
	}
	for(int c=0;c<paraCombNum;c++)
	{
		if(combWhetherUnique[rowNum][c]==true)
		{
			for(int i=0;i<t;i++)
			{
				if(index[cIndex[c][i]]==0)
					index[cIndex[c][i]]=1;
			}
		}
	}
	vector<int> newRow;
	for(int i=0;i<k;i++)
	{
		if(index[i]==1)
			newRow.insert(newRow.begin()+i,MCA[rowNum][i]);
		else
			newRow.insert(newRow.begin()+i,0);
	}
	changeFlag=false;
	int rowRT=getRuntime(MCA[rowNum]);
	curRT=rowRT;
	dfs2(newRow,index,0);
	if(changeFlag==true)
	{
		sumRT=sumRT-rowRT+curRT;
		vector<int> temp=MCA[rowNum];
		MCA[rowNum]=lessRTR;
		updateCIR(1,temp,rowNum,0,false);
		updateCIR(2,lessRTR,rowNum,0,false);
	}
	delete[] index;
}
void initGroup(int rowNum,int* index)
{
	group.clear();
	cMinRT.clear();
	NP=50;
	gNum=50;
	F=0.2;
	CR=0.2;
	for(int i=0;i<NP;i++)
	{
		vector<int> temp;
		for(int j=0;j<k;j++)
		{
			if(index[j]==1)
				temp.insert(temp.begin()+j,MCA[rowNum][j]);
			else
				temp.insert(temp.begin()+j,rand()%v[j]);
		}
		group.insert(group.begin()+i,temp);
		groupPro.insert(groupPro.begin()+i,getRuntime(temp));
	}
}
void findMin(int* index)
{
	int flag=0;
	for(int i=0;i<k;i++)
	{
		if(index[i]==0)
		{
			flag=1;
			break;
		}
	}
	if(flag==0)
	{
		cMinRTV=100000000;
		cMinRT=group[0];
		return;
	}
	int i=0;
	int minIndex=0;
	for(int j=1;j<NP;j++)
	{
		if(groupPro[j]<groupPro[minIndex])
		{
			minIndex=j;
		}
	}
	int minRuntime;
	int equalTime=0;
	while(i++<gNum)
	{
		minRuntime=groupPro[minIndex];
		for(int crow=0;crow<NP;crow++)
		{
			//mutation
			vector<int> candidateIn,fCandiIn;
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
			//float t=exp(1.0-(float)gNum/(float)(gNum+1-(i+1)));
			//float f=F*pow(2.0,t);
			for(int i=0;i<k;i++)
			{
				float point=(float)(group[r1][i]+F*(group[r2][i]-group[r3][i]))-(int)(group[r1][i]+F*(group[r2][i]-group[r3][i]));
				int temp=abs((int)(group[r1][i]+F*(group[r2][i]-group[r3][i])));
				int getLarger=rand()%10000;
				if((float)getLarger/10000.0<=point)
					temp=temp+1;
				temp=temp%v[i];
				candidateIn.insert(candidateIn.begin()+i,temp);
			}
			//crossover
			int jrand;
			do{jrand=rand()%k;}while(index[jrand]==1);
			for(int i=0;i<k;i++)
			{
				if((rand()%100<(CR*100)||i==jrand)&&index[i]!=1)
					fCandiIn.insert(fCandiIn.begin()+i,candidateIn[i]);
				else
					fCandiIn.insert(fCandiIn.begin()+i,group[crow][i]);
			}
			//selection
			//vector<int> fCandiIn=group[crow];
			int temp=getRuntime(fCandiIn);
			if(groupPro[crow]>temp)
			{
				group[crow]=fCandiIn;
				groupPro[crow]=temp;
				if(groupPro[crow]<groupPro[minIndex])
				{
					minIndex=crow;
				}
			}
		}
		if(minRuntime==groupPro[minIndex])
		{
			equalTime++;
			if(equalTime==10)
				break;
		}
		else
		{
			equalTime=0;
		}
	}
	cMinRT=group[minIndex];
	cMinRTV=groupPro[minIndex];
}
int calculateRowsRT(int rowNum,vector<int> fCandiIn,int tNum,int* ci,int d)
{
	int* index1=new int[k];
	for(int j=0;j<k;j++)
	{
		index1[j]=0;
	}
	for(int j=0;j<tNum;j++)
	{
		if(fCandiIn[j]==d)
		{
			for(int l=0;l<t;l++)
			{
				if(index1[cIndex[ci[j]][l]]==0)
					index1[cIndex[ci[j]][l]]=1;
			}
		}
	}
	initGroup(rowNum,index1);
	findMin(index1);
	return cMinRTV;
	delete[] index1;
}
int testValue;
int tNP1;
int tGNum1;
void findMinG(int rowNum,int degree,int curMinRT,int* index,int tNum,int* ci)//tNum is the unique combs' number
{
	if(degree==1)
	{
		initGroup(rowNum,index);
		findMin(index);
		candidateRows.insert(candidateRows.begin(),cMinRT);
		candidateRsRT=cMinRTV;
		return;
	}
	
	vector<vector<int>> group1;
	vector<int> groupPro1;
	int minIndex=0;
	int NP1=tNP1;
	int gNum1=tGNum1;
	float F1=0.2;
	float CR1=0.2;
	for(int i=0;i<NP1;i++)
	{
		vector<int> temp;
		for(int j=0;j<tNum;j++)
		{
			temp.insert(temp.begin()+j,rand()%degree);
		}
		int sumRT=0;
		for(int d=0;d<degree;d++)
		{
			calculateRowsRT(rowNum,temp,tNum,ci,d);
			sumRT=sumRT+cMinRTV;
		}
		group1.insert(group1.begin()+i,temp);
		groupPro1.insert(groupPro1.begin()+i,sumRT);
		if(sumRT<groupPro1[minIndex])
			minIndex=i;
	}
	int gCur=0;
	int minRuntime;
	int equalTime=0;
	while(gCur<gNum1)
	{
		gCur++;
		minRuntime=groupPro1[minIndex];
		for(int crow=0;crow<NP1;crow++)
		{
			//mutation
			vector<int> candidateIn,fCandiIn;
			int r1,r2,r3;
			do
			{
				r1=rand()%NP1;
			}while(r1==crow);
			do
			{
				r2=rand()%NP1;
			}while(r2==crow||r2==r1);
			do
			{
				r3=rand()%NP1;
			}while(r3==crow||r3==r1||r3==r2);
			//float t=exp(1.0-(float)gNum1/(float)(gNum1+1-(i+1)));
			//float f=F1*pow(2.0,t);
			for(int i=0;i<tNum;i++)
			{
				float point=(float)(group1[r1][i]+F1*(group1[r2][i]-group1[r3][i]))-(int)(group1[r1][i]+F1*(group1[r2][i]-group1[r3][i]));
				int temp=abs((int)(group1[r1][i]+F1*(group1[r2][i]-group1[r3][i])));
				int getLarger=rand()%10000;
				if((float)getLarger/10000.0<=point)
					temp=temp+1;
				temp=temp%degree;
				candidateIn.insert(candidateIn.begin()+i,temp);
			}
			//crossover
			int jrand;
			jrand=rand()%tNum;
			for(int i=0;i<tNum;i++)
			{
				if(rand()%100<(CR1*100)||i==jrand)
					fCandiIn.insert(fCandiIn.begin()+i,candidateIn[i]);
				else
					fCandiIn.insert(fCandiIn.begin()+i,group1[crow][i]);
			}
			//selection
			//vector<int> fCandiIn=group1[crow];
			int sumRT=0;
			for(int d=0;d<degree;d++)
			{
				calculateRowsRT(rowNum,fCandiIn,tNum,ci,d);
				sumRT=sumRT+cMinRTV;
				if(sumRT>groupPro1[crow])
					break;
			}
			if(groupPro1[crow]>sumRT)
			{
				group1[crow]=fCandiIn;
				groupPro1[crow]=sumRT;
				if(sumRT<groupPro1[minIndex])
					minIndex=crow;
			}
		}
		if(minRuntime==groupPro1[minIndex])
		{
			equalTime++;
			if(equalTime==10)
				break;
		}
		else
		{
			equalTime=0;
		}
	}
	if(groupPro1[minIndex]>=candidateRsRT)
		return;
	candidateRsRT=0;
	candidateRows.clear();
	for(int d=0;d<degree;d++)
	{
		int* index1=new int[k];
		for(int j=0;j<k;j++)
		{
			index1[j]=0;
		}
		for(int j=0;j<tNum;j++)
		{
			if(group1[minIndex][j]==d)
			{
				for(int l=0;l<t;l++)
				{
					if(index1[cIndex[ci[j]][l]]==0)
						index1[cIndex[ci[j]][l]]=1;
				}
			}
		}
		initGroup(rowNum,index1);
		findMin(index1);
		candidateRows.insert(candidateRows.begin()+d,cMinRT);
		candidateRsRT=candidateRsRT+cMinRTV;
		delete[] index1;
	}
}
void replaceOneRow2(int rowNum,int rows)
{
	candidateRows.clear();
	candidateRsRT=100000000;
	int* index=new int[k];
	for(int i=0;i<k;i++)
	{
		index[i]=0;
	}
	if(combNumInRow[rowNum]==0)
		return;
	int* ci=new int[combNumInRow[rowNum]];
	int tNum=0;
	for(int c=0;c<paraCombNum;c++)
	{
		if(combWhetherUnique[rowNum][c]==true)
		{
			for(int i=0;i<t;i++)
			{
				if(index[cIndex[c][i]]==0)
					index[cIndex[c][i]]=1;
			}
			ci[tNum]=c;
			tNum++;
		}
	}
	int rowUCNum=getRuntime(MCA[rowNum]);
	//for(int i=1;i<(tNum+1)/5;i++)
	for(int i=1;i<3;i++)
	{
		findMinG(rowNum,i,rowUCNum,index,tNum,ci);
	}
	delete[] index;
	delete[] ci;
	if(candidateRsRT>=rowUCNum)
		return;
	int size=candidateRows.size();
	sumRT=sumRT-rowUCNum+candidateRsRT;
	vector<int> temp=MCA[rowNum];
	MCA[rowNum]=candidateRows[0];
	updateCIR(1,temp,rowNum,0,false);
	updateCIR(2,candidateRows[0],rowNum,0,false);
	for(int i=1;i<size;i++)
	{
		addOneRow(N,2,candidateRows[i]);
		updateCIR(3,MCA[N],0,N,false);
		N=N+1;
	}
}
int canculateMinRT(vector<int> row,int index,int oldV)
{
	int indexV=oldV;
	for(int i=0;i<v[index];i++)
	{
		if(i==oldV)
			continue;
		row[index]=i;
		int temp=getRuntime(row);
		if(cMinRTV>temp)
		{
			cMinRTV=temp;
			indexV=i;
		}
	}
	return indexV;
}
bool findMinWithSA(vector<int> row,int* index,int RT)
{
	cMinRT=row;
	cMinRTV=RT;
	int flag=0;
	for(int i=0;i<k;i++)
	{
		if(index[i]==0)
		{
			flag=1;
			break;
		}
	}
	if(flag==0)
		return false;
	double T0=4.0,Tf=pow(10,-10),coolingFactor=0.99;
	int n=row.size();
	int V=n,frozenFactor=11;
	int L;
	T=T0;
	int count=0;
	cMinRT=row;
	cMinRTV=RT;
	int bestRT=RT;
	while(count<frozenFactor&&T>Tf)
	{
		//L=n*k*V*V;
		L=n*k*V;
		int storeCurRT=bestRT;
		for(int i=0;i<L;i++)
		{
			int j;
			do{j=rand()%k;}while(index[j]==1);
			int newV=canculateMinRT(row,j,row[j]);
			if(newV==row[j])
				continue;
			else
			{
				row[j]=newV;
				bestRT=cMinRTV;
			}
		}
		T=T*coolingFactor;
		if(bestRT==storeCurRT)
			count++;
		else
			count=0;
	}
	if(cMinRTV==RT)
	{
		return false;
	}
	else
	{
		cMinRT=row;
		return true;
	}
}
int calculateRowsRT2(int rowNum,vector<int> rows,int tNum,int* ci,int d)
{
	int* index1=new int[k];
	for(int j=0;j<k;j++)
	{
		index1[j]=0;
	}
	for(int j=0;j<tNum;j++)
	{
		if(rows[j]==d)
		{
			for(int l=0;l<t;l++)
			{
				if(index1[cIndex[ci[j]][l]]==0)
					index1[cIndex[ci[j]][l]]=1;
			}
		}
	}
	vector<int> temp;
	for(int j=0;j<k;j++)
	{
		if(index1[j]==1)
			temp.insert(temp.begin()+j,MCA[rowNum][j]);
		else
			temp.insert(temp.begin()+j,rand()%v[j]);
	}
	int RT=getRuntime(temp);
	findMinWithSA(temp,index1,RT);
	delete[] index1;
	return cMinRTV;
}
bool findMinWithSAG(int rowNum,int degree,int RT,int* index,int tNum,int* ci)
{
	if(degree==1)
	{
		if(findMinWithSA(MCA[rowNum],index,RT)==true)
		{
			candidateRows.insert(candidateRows.begin(),cMinRT);
			candidateRsRT=cMinRTV;
			return true;
		}
		else
		{
			candidateRows.insert(candidateRows.begin(),cMinRT);
			candidateRsRT=cMinRTV;
			return false;
		}
	}
	vector<int> rows;
	for(int j=0;j<tNum;j++)
	{
		rows.insert(rows.begin()+j,rand()%degree);
	}
	double T0=4.0,Tf=pow(10,-10),coolingFactor=0.15;
	int n=rows.size();
	int V=n,frozenFactor=11;
	int L;
	T=T0;
	int count=0;
	int sumRT=0;
	for(int d=0;d<degree;d++)
	{
		calculateRowsRT2(rowNum,rows,tNum,ci,d);
		sumRT=sumRT+cMinRTV;
	}
	int candidateRsRT1=sumRT;
	vector<int> cRows=rows;
	int bestRT=sumRT;
	while(count<frozenFactor&&T>Tf)
	{
		//L=n*k*V*V;
		L=5;
		int storeCurRT=bestRT;
		for(int i=0;i<L;i++)
		{
			int j=rand()%tNum;
			int oldV=rows[j];
			for(int m=0;m<degree;m++)
			{
				if(m==oldV)
					continue;
				rows[j]=m;
				int sumRT=0;
				for(int d=0;d<degree;d++)
				{
					calculateRowsRT2(rowNum,rows,tNum,ci,d);
					sumRT=sumRT+cMinRTV;
					if(sumRT>candidateRsRT1)
						break;
				}
				if(candidateRsRT1>sumRT)
				{
					candidateRsRT1=sumRT;
					cRows=rows;
					if(sumRT<bestRT)
						bestRT=sumRT;
				}
			}
		}
		T=T*coolingFactor;
		if(bestRT==storeCurRT)
			count++;
		else
			count=0;
	}
	if(candidateRsRT1>=candidateRsRT)
		return false;
	candidateRsRT=0;
	candidateRows.clear();
	for(int d=0;d<degree;d++)
	{
		int* index1=new int[k];
		for(int j=0;j<k;j++)
		{
			index1[j]=0;
		}
		for(int j=0;j<tNum;j++)
		{
			if(cRows[j]==d)
			{
				for(int l=0;l<t;l++)
				{
					if(index1[cIndex[ci[j]][l]]==0)
						index1[cIndex[ci[j]][l]]=1;
				}
			}
		}
		vector<int> temp;
		for(int j=0;j<k;j++)
		{
			if(index1[j]==1)
				temp.insert(temp.begin()+j,MCA[rowNum][j]);
			else
				temp.insert(temp.begin()+j,rand()%v[j]);
		}
		int RT=getRuntime(temp);
		if(findMinWithSA(temp,index1,RT)==true)
			candidateRows.insert(candidateRows.begin()+d,cMinRT);
		else
			candidateRows.insert(candidateRows.begin()+d,temp);
		candidateRsRT=candidateRsRT+cMinRTV;
	}
	return true;
}
void replaceOneRow3(int rowNum)
{
	candidateRows.clear();
	candidateRsRT=100000000;
	int* index=new int[k];
	for(int i=0;i<k;i++)
	{
		index[i]=0;
	}
	//if(combNumInRow[rowNum]==0)
		//return;
	int* ci=new int[combNumInRow[rowNum]];
	int tNum=0;
	for(int c=0;c<paraCombNum;c++)
	{
		if(combWhetherUnique[rowNum][c]==true)
		{
			for(int i=0;i<t;i++)
			{
				if(index[cIndex[c][i]]==0)
					index[cIndex[c][i]]=1;
			}
			ci[tNum]=c;
			tNum++;
		}
	}
	int rowUCNum=getRuntime(MCA[rowNum]);
	//findMinWithSAG(rowNum,1,rowUCNum,index,tNum,ci);
	for(int i=1;i<(tNum+1)/5;i++)
	{
		findMinWithSAG(rowNum,i,rowUCNum,index,tNum,ci);
		//cout<<i<<" "<<N<<" "<<candidateRows.size()<<endl;
		//verifyMCA();
	}
	if(candidateRsRT>=rowUCNum)
		return;
	int size=candidateRows.size();
	sumRT=sumRT-rowUCNum+candidateRsRT;
	vector<int> temp=MCA[rowNum];
	MCA[rowNum]=candidateRows[0];
	updateCIR(1,temp,rowNum,0,false);
	updateCIR(2,candidateRows[0],rowNum,0,false);
	for(int i=1;i<size;i++)
	{
		addOneRow(N,2,candidateRows[i]);
		updateCIR(3,MCA[N],0,N,false);
		N=N+1;
	}
}
void printMCA(ofstream &output)
{
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<k;j++)
		{
			//printf("%d ",MCA[i][j]);
			output<<MCA[i][j];
		}
		//printf("\n");
		output<<endl;
	}
	sumRT=0;
	for(int i=0;i<N;i++)
	{
		sumRT=sumRT+getRuntime(MCA[i]);
	}
	printf("totalCombNum is %d,N is %d,totalRuntime is %d.\n",totalCombNum,N,sumRT);
	output<<"totalCombNum,N,totalRuntime are "<<totalCombNum<<" "<<N<<" "<<sumRT<<endl<<endl<<endl;
	verifyMCA();

	if(deleteZRow()==true)
	{
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<k;j++)
			{
				//printf("%d ",MCA[i][j]);
				output<<MCA[i][j];
			}
			//printf("\n");
			output<<endl;
		}
		sumRT=0;
		for(int i=0;i<N;i++)
		{
			sumRT=sumRT+getRuntime(MCA[i]);
		}
		printf("totalCombNum is %d,N is %d,totalRuntime is %d.\n",totalCombNum,N,sumRT);
		output<<"totalCombNum,N,totalRuntime are "<<totalCombNum<<" "<<N<<" "<<sumRT<<endl<<endl<<endl;
		verifyMCA();
	}
	int storeN=N;
	vector<vector<int>> tMCA=MCA;
	vector<int> tCombNumInRow=combNumInRow;
	vector<int*> tCombIndex=combIndex;
	vector<vector<bool>> tCombWhetherUnique=combWhetherUnique;
	vector<int*> tCIndex=cIndex;
	vector<vector<key>>	tCombInRow=combInRow;
	unordered_map <key,int,hash_class,equal_class> tG=G;
	tNP1=15;tGNum1=15;
	int tempSum;
	//do
	//{
		tempSum=sumRT;
		for(int i=0;i<N;i++)
		{
			replaceOneRow2(i,2);
			//verifyMCA();
		}
	//}while(sumRT<tempSum);
	deleteZRow();
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<k;j++)
		{
			//printf("%d ",MCA[i][j]);
			output<<MCA[i][j];
		}
		//printf("\n");
		output<<endl;
	}
	sumRT=0;
	for(int i=0;i<N;i++)
	{
		sumRT=sumRT+getRuntime(MCA[i]);
	}
	printf("totalCombNum is %d,N is %d,totalRuntime is %d.\n",totalCombNum,N,sumRT);
	output<<"totalCombNum,N,totalRuntime are "<<totalCombNum<<" "<<N<<" "<<sumRT<<endl<<endl<<endl;
	verifyMCA();
	/*int count=0;
	while(count++<4)
	{
		tNP1=tNP1+5;
		//tGNum1=tGNum1+5;
	N=storeN;
	MCA=tMCA;
	G=tG;
	combNumInRow=tCombNumInRow;
	combIndex=tCombIndex;
	combWhetherUnique=tCombWhetherUnique;
	cIndex=tCIndex;
	combInRow=tCombInRow;
	//do
	//{
		tempSum=sumRT;
		for(int i=0;i<N;i++)
		{
			replaceOneRow2(i,2);
			//verifyMCA();
		}
	//}while(sumRT<tempSum);
	deleteZRow();
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<k;j++)
		{
			//printf("%d ",MCA[i][j]);
			output<<MCA[i][j];
		}
		//printf("\n");
		output<<endl;
	}
	sumRT=0;
	for(int i=0;i<N;i++)
	{
		sumRT=sumRT+getRuntime(MCA[i]);
	}
	printf("totalCombNum is %d,N is %d,totalRuntime is %d.\n",totalCombNum,N,sumRT);
	output<<"totalCombNum,N,totalRuntime are "<<totalCombNum<<" "<<N<<" "<<sumRT<<endl<<endl<<endl;
	verifyMCA();
	}
	/*N=storeN;
	MCA=tMCA;
	G=tG;
	combNumInRow=tCombNumInRow;
	combIndex=tCombIndex;
	combWhetherUnique=tCombWhetherUnique;
	cIndex=tCIndex;
	combInRow=tCombInRow;
	//do
	//{
		tempSum=sumRT;
		for(int i=0;i<N;i++)
		{
			replaceOneRow2(i,3);
			//verifyMCA();
		}
	//}while(sumRT<tempSum);
	deleteZRow();
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<k;j++)
		{
			//printf("%d ",MCA[i][j]);
			output<<MCA[i][j];
		}
		//printf("\n");
		output<<endl;
	}
	sumRT=0;
	for(int i=0;i<N;i++)
	{
		sumRT=sumRT+getRuntime(MCA[i]);
	}
	printf("totalCombNum is %d,N is %d,totalRuntime is %d.\n",totalCombNum,N,sumRT);
	output<<"totalCombNum,N,totalRuntime are "<<totalCombNum<<" "<<N<<" "<<sumRT<<endl<<endl<<endl;
	verifyMCA();*/

	/*if(deleteZRow()==true)
	{
		for(int i=0;i<N;i++)
		{
			for(int j=0;j<k;j++)
			{
				//printf("%d ",MCA[i][j]);
				output<<MCA[i][j];
			}
			//printf("\n");
			output<<endl;
		}
		sumRT=0;
		for(int i=0;i<N;i++)
		{
			sumRT=sumRT+getRuntime(MCA[i]);
		}
		printf("totalCombNum is %d,N is %d,totalRuntime is %d.\n",totalCombNum,N,sumRT);
		output<<"totalCombNum,N,totalRuntime are "<<totalCombNum<<" "<<N<<" "<<sumRT<<endl<<endl<<endl;
		verifyMCA();
	}*/
}
int main()
{
	ifstream input;
	input.open("testcase.txt");
	while(input.peek()!=EOF)
	{
		v.clear();

		ofstream output;
		string name;
		stringstream ss;
		printf("\n\nNew\n");
		output<<endl<<endl<<"New"<<endl;
		printf("Input parameter t k p v1 v2 ... vk:\n");
		output<<"Input parameter t k p v1 v2 ... vk:"<<endl;
		//scanf("%d%d%d",&t,&k,&p);
		input>>t>>k>>p;
		ss<<t<<0<<k<<0<<p<<0;
		if(k==1)
		{
			printf("error!\n");
			return 0;
		}
		int* vt=new int[k];
		for(int i=0;i<k;i++)
		{
			int temp;
			//scanf("%d",&temp);
			input>>temp;
			ss<<temp<<0;
			if(temp==0)
			{
				printf("error!\n");
				return 0;
			}
			v.insert(v.begin()+i,temp);
			vt[i]=temp;
		}
		name=ss.str();
		name=name.append(".txt");
		output.open(name);
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
		int storeN=N;
		int countCS=0;
		delete[] vt;
		while(countCS<3)
		{
			N=storeN;
			countCS++;
			//initialize parameters and MCA
			totalCombNum=0;
			curCombNum=0;
			paraCombNum=0;
			findMCA=false;
			MCA.clear();
			combinations.clear();
			for(int i=0;i<combIndex.size();i++)
			{
				delete[] combIndex[i];
			}
			combIndex.clear();
			combNumInRow.clear();
			combWhetherUnique.clear();
			
			runtime.clear();
			for(int i=0;i<cIndex.size();i++)
			{
				delete[] cIndex[i];
			}
			cIndex.clear();
			lessRTR.clear();
			group.clear();
			groupPro.clear();
			combInRow.clear();
			G.clear();
			
			initMCA();
			genCombs();
			getAllCaseRuntime(2);
			for(int i=0;i<N;i++)
			{
				addCombToG(MCA[i]);
			}
			if(fitness(1,0,MCA[0])==0)
			{
				findMCA=true;
				printMCA(output);
				continue;
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
							findMCA=true;
							printMCA(output);
							break;
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
				if(findMCA==true)
					break;
				addOneRow(N,1,MCA[0]);
				updateCIR(3,MCA[N],0,N,false);
				N=N+1;
				if(fitness(1,0,MCA[0])==0)
				{
					findMCA=true;
					printMCA(output);
					break;
				}
			}
		}
		output.close();
	}
	input.close();
	return 0;
}