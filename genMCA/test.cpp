#include<iostream>
#include<fstream>
#include<vector>
#include<string.h>
#include "t.h"
using namespace std;
int f(int a[2])
{
	a[0] = 1;
	return a[0];
}
int main()
{
	//jenny
	/*char num[2];
	char ch;
	vector<int> row;
	int rowNum=0;
	int sum=0;
	int ab;
	ifstream input;
	input.open("result.txt");
	ofstream output;
	output.open("data.txt");
	while(input.peek()!=EOF)
	{
		input.get(ch);
		if(ch=='\n')
		{
			cout<<rowNum<<"  "<<sum<<endl;
			output<<rowNum<<"  "<<sum<<endl;
			sum=0;
			rowNum=0;
			continue;
		}
		int i=0;
		rowNum++;
		row.clear();
		while(ch!='\n')
		{
			input>>ab;
			input>>num;
			row.insert(row.begin()+i,num[0]-'a');
			i++;
			input.get(ch);
		}
		for(int j=0;j<row.size()-1;j++)
		{
			sum=sum+row[j]*row[j+1];
		}
	}
	input.close();*/

	//allpairs
	/*char num[5];
	char ch;
	vector<int> row;
	int rowNum=0;
	int sum=0;
	ifstream input;
	input.open("result.txt");
	ofstream output;
	output.open("data.txt");
	while(input.peek()!=EOF)
	{
		input.get(ch);
		if(ch=='\n')
		{
			cout<<rowNum<<"  "<<sum<<endl;
			output<<rowNum<<"  "<<sum<<endl;
			sum=0;
			rowNum=0;
			continue;
		}
		int i=0;
		rowNum++;
		row.clear();
		char temp[3];
		if(rowNum>=10)
			input>>temp;
		while(ch!='\n')
		{
			input>>num;
			input.get(ch);
			if(ch=='\n')
				break;
			if(num[0]!='~')
				row.insert(row.begin()+i,num[0]-'0');
			else
				row.insert(row.begin()+i,num[1]-'0');
			i++;
			
		}
		for(int j=0;j<row.size()-1;j++)
		{
			sum=sum+row[j]*row[j+1];
		}
	}
	input.close();*/

	//acts
	char ch;
	vector<int> row;
	int rowNum=0;
	int sum=0;
	int value;
	ifstream input;
	input.open("result.txt");
	ofstream output;
	output.open("data.txt");
	while(input.peek()!=EOF)
	{
		input.get(ch);
		if(ch=='\n')
		{
			cout<<rowNum<<"  "<<sum<<endl;
			output<<rowNum<<"  "<<sum<<endl;
			sum=0;
			rowNum=0;
			continue;
		}
		rowNum++;
		row.clear();
		row.insert(row.begin(),ch-'0');
		int i=1;
		while(ch!='\n')
		{
			input.get(ch);
			if(ch=='\n')
				break;
			input.get(ch);
			row.insert(row.begin()+i,ch-'0');
			i++;
		}
		double temp = 6.7713-5.5377*row[3]+(1-row[3])*(840.68*row[4]-835.51*row[4]*row[7]-102.27*(1-row[7])*(row[0]-1)*row[0]*row[4]);
		sum=sum+temp;
		/*for(int j=0;j<row.size()-1;j++)
		{
			sum=sum+row[j]*row[j+1];
		}*/
	}
	cout<<sum<<endl;
	input.close();

	//pict
	/*char ch;
	vector<int> row;
	int rowNum=0;
	int sum=0;
	int value;
	ifstream input;
	input.open("result.txt");
	ofstream output;
	output.open("data.txt");
	while(input.peek()!=EOF)
	{
		input.get(ch);
		if(ch=='\n')
		{
			cout<<rowNum<<"  "<<sum<<endl;
			output<<rowNum<<"  "<<sum<<endl;
			sum=0;
			rowNum=0;
			continue;
		}
		rowNum++;
		row.clear();
		row.insert(row.begin(),ch-'0');
		int i=1;
		while(ch!='\n')
		{
			input.get(ch);
			if(ch=='\n')
				break;
			input.get(ch);
			row.insert(row.begin()+i,ch-'0');
			i++;
		}
		for(int j=0;j<row.size()-1;j++)
		{
			sum=sum+row[j]*row[j+1];
		}
	}
	input.close();*/

	//testcover
	/*char num[10];
	char ch;
	vector<int> row;
	int rowNum=0;
	int sum=0;
	ifstream input;
	input.open("result.txt");
	ofstream output;
	output.open("data.txt");
	while(input.peek()!=EOF)
	{
		input.get(ch);
		if(ch=='\n')
		{
			cout<<rowNum<<"  "<<sum<<endl;
			output<<rowNum<<"  "<<sum<<endl;
			sum=0;
			rowNum=0;
			continue;
		}
		int i=0;
		rowNum++;
		row.clear();
		char temp[10];
		if(rowNum>=10)
			input>>temp;
		while(ch!='\n')
		{
			input>>num;
			input.get(ch);
			if(ch=='\n')
				break;
			row.insert(row.begin()+i,num[0]-'0');
			i++;
			
		}
		for(int j=0;j<row.size()-1;j++)
		{
			sum=sum+row[j]*row[j+1];
		}
	}
	input.close();*/
	//protest
	/*char ch;
	vector<int> row;
	int rowNum=0;
	int sum=0;
	int value;
	ifstream input;
	input.open("result.txt");
	ofstream output;
	output.open("data.txt");
	while(input.peek()!=EOF)
	{
		input.get(ch);
		if(ch=='\n')
		{
			cout<<rowNum<<"  "<<sum<<endl;
			output<<rowNum<<"  "<<sum<<endl;
			sum=0;
			rowNum=0;
			continue;
		}
		rowNum++;
		row.clear();
		row.insert(row.begin(),ch-'0');
		int i=1;
		while(ch!='\n')
		{
			input.get(ch);
			if(ch=='\n')
				break;
			input.get(ch);
			row.insert(row.begin()+i,ch-'0');
			i++;
		}
		for(int j=0;j<row.size()-1;j++)
		{
			sum=sum+row[j]*row[j+1];
		}
	}
	input.close();*/
	/*ofstream output("data.txt");
	int a;
	for(a=1001;a<4001;a++)
		output<<a<<endl;*/
}