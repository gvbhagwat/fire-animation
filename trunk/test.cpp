#include<iostream>
#include<climits>
using namespace std;


inline unsigned int randhash(unsigned int seed)
{
	unsigned int i = (seed^0xA3c59AC3u)*2654435769u;
	i^=(i>>16);
	i*=2654435769u;
	i^=(i>>16);
	i*=2654435769u;
	return i;

}	

inline double randhashd(unsigned int seed, double a, double b)
{
	return (b-a)*randhash(seed)/(double)UINT_MAX + a;
}

int main()
{
	for(int i = 1; i < 41; i++)
		for(int j = 1; j < 30; j++)
		{
			double x = randhashd(i*j^2, 0 , 1)+i;
			double y = randhashd(i^2*j+1, 0, 1)+j;
			cout<<"x = "<<x<<" y = "<<y<<endl;
		}
	
	return 0;
}		
